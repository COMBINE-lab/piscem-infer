//! Locating and validating a bulk RAD input, for both the legacy (v1) prelude
//! and the versioned (v2, `RAD_FILE` magic + tag roles) prelude.

use std::path::{Path, PathBuf};

/// Resolve the `-i` argument to the RAD file and its `map_info.json`.
///
/// `piscem map-bulk -o <stem>` writes `<stem>.rad` and `<stem>.map_info.json`
/// by appending to the stem's file name, so do the same here; replacing the
/// extension instead would turn the stem `sample.rep1` into `sample.rad`. The
/// path of the `.rad` file itself is accepted as well.
pub(crate) fn input_paths(input: &Path) -> (PathBuf, PathBuf) {
    let stem = if input.extension().is_some_and(|e| e == "rad") {
        input.with_extension("")
    } else {
        input.to_path_buf()
    };
    (
        append_to_file_name(&stem, ".rad"),
        append_to_file_name(&stem, ".map_info.json"),
    )
}

/// Append `suffix` to the last component of `stem` (as `PathBuf::add_extension`
/// does), so that e.g. `out/` becomes `out.rad` rather than `out/.rad`.
fn append_to_file_name(stem: &Path, suffix: &str) -> PathBuf {
    let mut name = stem.file_name().unwrap_or_default().to_os_string();
    name.push(suffix);
    stem.with_file_name(name)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn input_paths_append_to_the_stem() {
        let p = |s: &str| PathBuf::from(s);
        let (rad, info) = input_paths(Path::new("out/SRR1.rep1"));
        assert_eq!(rad, p("out/SRR1.rep1.rad"));
        assert_eq!(info, p("out/SRR1.rep1.map_info.json"));
        // the .rad file itself resolves the same way (and never to e.g. a
        // sibling out/SRR1.rad)
        assert_eq!(input_paths(Path::new("out/SRR1.rep1.rad")), (rad, info));
        assert_eq!(
            input_paths(Path::new("sample")),
            (p("sample.rad"), p("sample.map_info.json"))
        );
        assert_eq!(
            input_paths(Path::new("sample.rad")),
            (p("sample.rad"), p("sample.map_info.json"))
        );
        assert_eq!(
            input_paths(Path::new("bulk_0.23.0")),
            (p("bulk_0.23.0.rad"), p("bulk_0.23.0.map_info.json"))
        );
        // dots in a parent directory don't matter
        assert_eq!(
            input_paths(Path::new("run.1/sample")),
            (p("run.1/sample.rad"), p("run.1/sample.map_info.json"))
        );
        // like `piscem map-bulk -o out/`, which writes out.rad
        assert_eq!(
            input_paths(Path::new("out/")),
            (p("out.rad"), p("out.map_info.json"))
        );
    }
}
