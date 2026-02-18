# Feature Tracking

## Automatic Library Type Detection

**Status:** Implemented (not yet released)

**Summary:** Users can specify `--lib-type auto` to automatically detect the library type from a sample of mapped reads, instead of specifying it explicitly (SF, ISF, SR, ISR, U, IU).

**How it works:**
- Samples up to 10,000 mapped reads (configurable via `--auto-detect-samples`)
- Counts orientation frequencies (Forward, Reverse, FR, RF, etc.)
- Uses salmon-style heuristics (30%/70% thresholds) to classify as stranded or unstranded
- Paired-end: ISF / ISR / IU; Single-end: SF / SR / U
- Warns if ≥5% wrong-strand mappings (stranded) or >5% strand bias (unstranded)
- Detected type recorded in `.meta_info.json` as `inferred_lib_type`

**Files modified:**
- `src/utils/map_record_types.rs` — `OrientationCounts`, `detect_library_type()`, `check_strand_warnings()`, `Display` for `LibraryType`
- `src/prog_opts.rs` — `LibTypeArg` enum, `auto_detect_samples` CLI arg
- `src/process_rad.rs` — `detect_lib_type_from_sample()`, resolution in `process_bulk_dispatch()`
