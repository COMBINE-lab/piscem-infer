use clap::Parser;

use tracing::Level;

mod fld;
mod process_rad;
mod prog_opts;
mod utils;
use prog_opts::{Cli, Commands};
use utils::eq_maps::EqMapType;

use crate::process_rad::process_bulk;

fn main() -> anyhow::Result<()> {
    let cli_args = Cli::parse();
    if cli_args.quiet {
        tracing_subscriber::fmt().with_max_level(Level::WARN).init();
    } else {
        tracing_subscriber::fmt().with_max_level(Level::INFO).init();
    }

    match cli_args.command {
        Commands::Quant(quant_opts) => {
            // if there is only 1 probability bin, then all fragments
            // fall in the same probability bin and we are using the basic
            // equivalence class definition. Otherwise, we are using the
            // range factorized equivalence class definition.
            let eqc_type = if quant_opts.factorized_eqc_bins > 1 {
                EqMapType::RangeFactorizedEqMap
            } else {
                EqMapType::BasicEqMap
            };
            // set the
            crate::utils::eq_maps::NUM_BINS
                .set(quant_opts.factorized_eqc_bins as f64)
                .expect("NUM_BINS should not yet have been initalized");
            process_bulk(quant_opts, eqc_type)?
        }
    }
    Ok(())
}
