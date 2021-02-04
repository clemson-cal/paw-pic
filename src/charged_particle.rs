use crate::four_momentum::FourMomentum;
use crate::three_vector::ThreeVector;

/**
 * Describes a relativistic charged particle
 */
#[derive(Clone)]
struct ChargedParticle {
    pub charge: f64,
    pub position: ThreeVector,
    pub momentum: FourMomentum,
}

// ============================================================================
impl ChargedParticle {
    /**
     * Move the particle forward in time using the classic "Boris Push" method.
     */
    pub fn boris_push(&self, electric: ThreeVector, magnetic: ThreeVector, dt: f64) -> Self {
        let e = self.charge;
        let m = self.momentum.rest_mass();
        let h = e / m * dt;

        let u_nmh = self.momentum.gamma_beta_vector();
        let u_minus = u_nmh + electric * (0.5 * h);
        let gamma_n = (1.0 + u_minus.squared()).sqrt();
        let t = magnetic * (0.5 * h / gamma_n);
        let s = t / (1.0 + t.squared()) * 2.0;
        let u_prime = u_minus + u_minus.cross(&t);
        let u_plus = u_minus + u_prime.cross(&s);
        let u_nph = u_plus + electric * (0.5 * h);
        let gamma_nph = (1.0 + u_nph.squared()).sqrt();
        let new_pos = self.position + u_nph * dt / gamma_nph;
        let new_mom = FourMomentum(gamma_nph, u_nph.0, u_nph.1, u_nph.2) * m;

        ChargedParticle {
            charge: e,
            position: new_pos,
            momentum: new_mom,
        }
    }

    /*
     * Move the particle forward in time using rk4 push method.
     */
    pub fn rk4_push(&self, electric: ThreeVector, magnetic: ThreeVector, dt: f64) -> Self {
        let e = self.charge;
        let m = self.momentum.rest_mass();
        let h = e / m * dt;
        let u0 = self.momentum.gamma_beta_vector();

        let dxdt_0 = self.momentum.velocity_vector();
        let dudt_0 = (electric + u0.cross(&magnetic) / self.momentum.lorentz_factor()) * (e / m);

        let x_one = self.position + dxdt_0 * dt * 0.5;
        let u_one = self.position + dudt_0 * dt * 0.5;
        let p_one = FourMomentum(self.momentum.0, u_one.0, u_one.1, u_one.2);

        let dxdt_1 = u_one / p_one.lorentz_factor();
        let dudt_1 = (electric + u_one.cross(&magnetic) / p_one.lorentz_factor()) * (e / m);

        let x_two = self.position + dxdt_1 * dt * 0.5;
        let u_two = self.position + dudt_1 * dt * 0.5;
        let p_two = FourMomentum(self.momentum.0, u_two.0, u_two.1, u_two.2);

        let dxdt_2 = u_two / p_two.lorentz_factor();
        let dudt_2 = (electric + u_two.cross(&magnetic) / p_two.lorentz_factor()) * (e / m);

        let x_three = self.position + dxdt_2 * dt;
        let u_three = self.position + dudt_2 * dt;
        let p_three = FourMomentum(self.momentum.0, u_three.0, u_three.1, u_three.2);

        let dxdt_3 = u_three / p_three.lorentz_factor();
        let dudt_3 = (electric + u_three.cross(&magnetic) / p_three.lorentz_factor()) * (e / m);

        let new_pos = self.position + (dxdt_0 + (dxdt_1 + dxdt_2) * 2.0 + dxdt_3) * (dt / 6.0);
        let new_u = u0 + (dudt_0 + (dudt_1 + dudt_2) * 2.0 + dudt_3) * (dt / 6.0);
        let new_mom = FourMomentum(self.momentum.0, new_u.0, new_u.1, new_u.2);

        ChargedParticle {
            charge: e,
            position: new_pos,
            momentum: new_mom,
        }
    }
}
