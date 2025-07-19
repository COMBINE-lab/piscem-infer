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
        (1.0 - self.cum_probs.get(i).unwrap_or(&1.0)) + f64::MIN_POSITIVE
    }
    fn cdf(&self, i: usize) -> f64 {
        self.cum_probs.get(i).unwrap_or(&1.0) + f64::MIN_POSITIVE
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
        self.cum_probs.get(i).unwrap_or(&1.0) + f64::MIN_POSITIVE
    }
}
