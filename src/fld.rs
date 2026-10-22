pub enum Fld {
    Empirical(EmpiricalFLD),
    Parametric(ParametricFLD),
}

pub trait FldPDF {
    fn pdf(&self, i: usize) -> f64;
    fn cdf(&self, i: usize) -> f64;
    #[allow(dead_code)]
    fn survival(&self, i: usize) -> f64;
}

pub struct ParametricFLD {
    mu: f64,
    sigma: f64,
    inv_sigma: f64,
    inv_denom: f64,
    cum_probs: Vec<f64>,
}

impl ParametricFLD {
    pub fn new(mu: f64, sigma: f64, upper: usize) -> Self {
        let inv_sigma = 1.0 / sigma;
        let denom_b = distrs::Normal::cdf(upper as f64, mu, sigma);
        let denom_a = distrs::Normal::cdf(0.0_f64, mu, sigma);
        let denom = denom_b - denom_a;
        let inv_denom = 1.0_f64 / denom;

        let cum_probs: Vec<f64> = (0..upper)
            .scan(0.0, |sum, x| {
                let px = inv_sigma * (distrs::Normal::pdf(x as f64, mu, sigma) * inv_denom);
                *sum += px;
                Some(*sum)
            })
            .collect();

        Self {
            mu,
            sigma,
            inv_sigma,
            inv_denom,
            cum_probs,
        }
    }
}

impl FldPDF for ParametricFLD {
    fn pdf(&self, i: usize) -> f64 {
        let x = i as f64;
        self.inv_sigma * (distrs::Normal::pdf(x, self.mu, self.sigma) * self.inv_denom)
    }
    fn survival(&self, i: usize) -> f64 {
        // same past-the-table fallback as `cdf` (see there)
        let last = self.cum_probs.last().unwrap_or(&1.0);
        (1.0 - self.cum_probs.get(i).unwrap_or(last)) + f64::MIN_POSITIVE
    }
    fn cdf(&self, i: usize) -> f64 {
        // `cum_probs` tops out at ~1/sigma rather than 1 (each mass carries an
        // extra 1/sigma, which cancels since these values are only compared
        // within a fragment), so lengths past the end of the table take its last
        // entry; falling back to 1.0 would weight them ~sigma-fold above any
        // in-range length.
        let last = self.cum_probs.last().unwrap_or(&1.0);
        self.cum_probs.get(i).unwrap_or(last) + f64::MIN_POSITIVE
    }
}

pub struct EmpiricalFLD {
    probs: Vec<f64>,
    cum_probs: Vec<f64>,
    epsilon: f64,
}

impl EmpiricalFLD {
    pub fn new(counts: Vec<u32>, eps: f64) -> Self {
        let tot_eps = eps * counts.len() as f64;
        let tot_mass = tot_eps + counts.iter().fold(0.0, |acc, x| acc + *x as f64);
        let probs: Vec<f64> = counts
            .iter()
            .map(|x| ((*x as f64) + eps) / tot_mass)
            .collect();
        let cum_probs = probs
            .iter()
            .scan(0.0, |sum, i| {
                *sum += i;
                Some(*sum)
            })
            .collect();
        Self {
            probs,
            cum_probs,
            epsilon: eps,
        }
    }
}

impl FldPDF for EmpiricalFLD {
    fn pdf(&self, i: usize) -> f64 {
        *self.probs.get(i).unwrap_or(&self.epsilon)
    }
    fn survival(&self, i: usize) -> f64 {
        (1.0 - self.cum_probs.get(i).unwrap_or(&1.0)) + self.epsilon
    }
    fn cdf(&self, i: usize) -> f64 {
        self.cum_probs.get(i).unwrap_or(&1.0) + self.epsilon
    }
}

#[cfg(test)]
mod tests {
    use super::{FldPDF, ParametricFLD};
    use crate::utils::em::{adjust_ref_lengths, conditional_means, conditional_means_from_params};

    // the upper bound used for the FLD and conditional means in `process_rad`
    const UPPER: usize = 65_536;

    #[test]
    fn parametric_cdf_is_flat_past_upper() {
        for (mu, sigma) in [(200.0, 80.0), (250.0, 25.0)] {
            let fld = ParametricFLD::new(mu, sigma, UPPER);
            // well past the mean the cdf has saturated; lengths at or past the end
            // of the table must get that same weight, not a sigma-fold larger one
            let saturated = fld.cdf(5_000);
            assert_eq!(fld.cdf(UPPER - 1), saturated);
            assert_eq!(fld.cdf(UPPER), saturated);
            assert_eq!(fld.cdf(100_000), saturated);
            assert!(fld.cdf(100) < fld.cdf(200) && fld.cdf(200) < saturated);
            assert_eq!(fld.survival(UPPER), fld.survival(UPPER - 1));
        }
    }

    #[test]
    fn adjust_ref_lengths_at_fld_upper_bound() {
        // with a flat FLD over 0..UPPER, cond_means[i] = i / 2
        let cond_means = conditional_means(&vec![1u32; UPPER]);
        assert_eq!(cond_means.len(), UPPER);
        let eff_lens = adjust_ref_lengths(&[65_535, 65_536, 65_537], &cond_means);
        assert_eq!(eff_lens, vec![32_767.5, 32_768.5, 32_769.5]);

        let cond_means = conditional_means_from_params(200.0, 80.0, UPPER);
        let tmean = *cond_means.last().unwrap();
        let eff_lens = adjust_ref_lengths(&[65_535, 65_536, 65_537], &cond_means);
        assert_eq!(
            eff_lens,
            vec![65_535.0 - tmean, 65_536.0 - tmean, 65_537.0 - tmean]
        );
    }
}
