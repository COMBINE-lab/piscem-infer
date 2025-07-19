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
            let eqc_type = if quant_opts.factorized_eqc_bins > 1 {
                EqMapType::RangeFactorizedEqMap
            } else {
                EqMapType::BasicEqMap
            };
            process_bulk(quant_opts, eqc_type)?
        }
    }
    Ok(())
}
