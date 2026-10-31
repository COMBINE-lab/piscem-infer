#!/usr/bin/env bash
set -euo pipefail

# Bump piscem-infer's version, record the release in CHANGELOG.md, commit, and
# (with --publish) publish to crates.io; then tag and push. Pushing the tag
# starts the release workflow (cargo-dist), which builds the binaries and
# creates the GitHub release. Run it on an up-to-date, clean `main`.

die() {
    echo "error: $*" >&2
    exit 1
}

usage() {
    cat <<'EOF'
Usage:
  ./bump_and_publish.sh <version> [--publish] [--dry-run] [--no-changelog]

Options:
  --publish       Publish to crates.io after bumping and committing
  --dry-run       Show what would be done without modifying files, creating
                  commits or tags, pushing, or publishing (the packaging check
                  still runs)
  --no-changelog  Do not add a CHANGELOG.md section (requires git-cliff otherwise)
  -h, --help      Show this help message

Without --publish the bump is still committed, tagged and pushed; only the
crates.io upload is skipped.

The new CHANGELOG.md section is generated from conventional commits by
git-cliff (see cliff.toml) and prepended to the file. Preview it, without
writing anything:

  git-cliff --tag v<version> --unreleased
EOF
}

print_cmd() {
    printf '+'
    printf ' %q' "$@"
    printf '\n'
}

run() {
    print_cmd "$@"
    if [[ "$DRY_RUN" == true ]]; then
        return 0
    fi
    "$@"
}

VERSION=""
PUBLISH=false
DRY_RUN=false
CHANGELOG_ENABLED=true

while [[ $# -gt 0 ]]; do
    case "$1" in
        --publish) PUBLISH=true ;;
        --dry-run) DRY_RUN=true ;;
        --no-changelog) CHANGELOG_ENABLED=false ;;
        -h|--help) usage; exit 0 ;;
        -*) die "unknown option: $1" ;;
        *)
            [[ -z "$VERSION" ]] || die "version specified more than once"
            VERSION="$1"
            ;;
    esac
    shift
done

[[ -n "$VERSION" ]] || { usage; exit 1; }

if ! [[ "$VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+([+-][0-9A-Za-z.-]+)*$ ]]; then
    die "version must look like X.Y.Z, optionally with prerelease/build suffixes"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

CARGO="Cargo.toml"
LOCKFILE="Cargo.lock"
CHANGELOG="CHANGELOG.md"
CLIFF_CONFIG="cliff.toml"
CRATE="piscem-infer"
TAG="v${VERSION}"

[[ -f "$CARGO" ]] || die "not found: $CARGO"
[[ -f "$LOCKFILE" ]] || die "not found: $LOCKFILE"
if [[ "$CHANGELOG_ENABLED" == true ]]; then
    [[ -f "$CLIFF_CONFIG" ]] || die "not found: $CLIFF_CONFIG (pass --no-changelog to skip the changelog)"
    [[ -f "$CHANGELOG" ]] || die "not found: $CHANGELOG"
    command -v git-cliff >/dev/null 2>&1 || die "git-cliff is not installed; install it (cargo binstall git-cliff) or pass --no-changelog"
fi

# A [patch] section would make the published crate build against something
# other than what was tested here; crates.io ignores it, so refuse instead.
if grep -q '^\[patch' "$CARGO"; then
    die "$CARGO has a [patch] section; depend on published crates before releasing"
fi

CURRENT_VERSION="$(sed -n 's/^version = "\(.*\)"/\1/p' "$CARGO" | head -1)"
[[ -n "$CURRENT_VERSION" ]] || die "could not determine the current version from $CARGO"
[[ "$CURRENT_VERSION" != "$VERSION" ]] || die "$CRATE is already at $VERSION"

if git rev-parse "$TAG" >/dev/null 2>&1; then
    die "tag $TAG already exists"
fi
if [[ -n "$(git status --porcelain)" ]]; then
    die "working tree is not clean; commit or stash existing changes first"
fi
git remote get-url origin >/dev/null 2>&1 || die "git remote 'origin' is not configured"

echo "Current version : $CURRENT_VERSION"
echo "New version     : $VERSION"
echo "Tag             : $TAG"
echo "Publish         : $([[ "$PUBLISH" == true ]] && echo yes || echo no)"
echo "Changelog       : $([[ "$CHANGELOG_ENABLED" == true ]] && echo "prepend a $TAG section" || echo no)"
echo "Dry-run         : $([[ "$DRY_RUN" == true ]] && echo yes || echo no)"
echo

if [[ "$DRY_RUN" == false ]]; then
    # Only the [package] version (the first `version = ` line), and the
    # crate's own entry in the lockfile.
    sed -i.bak "1,/^version = /s/^version = \".*\"/version = \"${VERSION}\"/" "$CARGO"
    rm -f "${CARGO}.bak"
    sed -i.bak "/^name = \"${CRATE}\"$/,/^dependencies = \\[$/s/^version = \".*\"/version = \"${VERSION}\"/" "$LOCKFILE"
    rm -f "${LOCKFILE}.bak"
    [[ "$(sed -n 's/^version = "\(.*\)"/\1/p' "$CARGO" | head -1)" == "$VERSION" ]] || die "version update failed"
else
    echo "Dry-run: would set the version to $VERSION in $CARGO and $LOCKFILE"
fi

# Check exactly what crates.io will package (and that it builds) before
# anything irreversible. --locked: the bump must not change any dependency.
cargo publish --dry-run --locked --allow-dirty

if [[ "$CHANGELOG_ENABLED" == true ]]; then
    # `--tag` labels the not-yet-tagged commits with the version being cut;
    # without it they would land under "Unreleased".
    echo "Prepending the $TAG section to $CHANGELOG"
    run git-cliff --tag "$TAG" --unreleased --prepend "$CHANGELOG"
fi

run git add "$CARGO" "$LOCKFILE"
if [[ "$CHANGELOG_ENABLED" == true ]]; then
    run git add "$CHANGELOG"
fi
run git commit -m "chore(release): bump ${CRATE} to v${VERSION}"

if [[ "$PUBLISH" == true ]]; then
    run cargo publish --locked
fi

run git tag -a "$TAG" -m "Release ${VERSION}"
run git push origin HEAD
run git push origin "$TAG"

echo
if [[ "$DRY_RUN" == true ]]; then
    echo "Dry-run complete"
else
    echo "Release complete for v${VERSION}"
fi
