use crate::four_momentum::FourMomentum;
use crate::three_vector::ThreeVector;

/**
 * Describes a relativistic charged particle
 */
#[derive(Debug, Clone)]
pub struct ChargedParticle {
    pub charge: f64,
    pub position: ThreeVector,
    pub momentum: FourMomentum,
}

// ============================================================================
impl ChargedParticle {
    /**
     * Move the particle forward in time using the classic "Boris Push" method.
     */
    pub fn total_energy(&self) -> f64 {
        let p = &self.momentum;
        p.rest_mass() * p.lorentz_factor()
    }
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
    pub fn rk4_push<F>(&self, field: F, time: f64, dt: f64) -> Self
    where
        F: Fn(ThreeVector, f64) -> (ThreeVector, ThreeVector),
    {
        let e = self.charge;
        let m = self.momentum.rest_mass();

        let deltas = |x: ThreeVector, u: ThreeVector, t: f64| -> (ThreeVector, ThreeVector, f64) {
            let (electric, magnetic) = field(x, t);
            let v = u / (1.0 + u.squared()).sqrt();
            let dxdt = v;
            let dudt = (electric + v.cross(&magnetic)) * (e / m);
            (dxdt * dt, dudt * dt, dt)
        };

        let t0 = time;
        let x0 = self.position;
        let u0 = self.momentum.gamma_beta_vector();

        let (dx0, du0, dt0) = deltas(x0, u0, t0);
        let x1 = x0 + dx0 * 0.5;
        let u1 = u0 + du0 * 0.5;
        let t1 = t0 + dt0 * 0.5;

        let (dx1, du1, dt1) = deltas(x1, u1, t1);
        let x2 = x1 + dx1 * 0.5;
        let u2 = u1 + du1 * 0.5;
        let t2 = t1 + dt1 * 0.5;

        let (dx2, du2, dt2) = deltas(x2, u2, t2);
        let x3 = x2 + dx2;
        let u3 = u2 + du2;
        let t3 = t2 + dt2;

        let (dx3, du3, dt3) = deltas(x3, u3, t3);

        let x = x0 + (dx0 + dx1 * 2.0 + dx2 * 2.0 + dx3) / 6.0;
        let u = u0 + (du0 + du1 * 2.0 + du2 * 2.0 + du3) / 6.0;
        let t = t0 + (dt0 + dt1 * 2.0 + dt2 * 2.0 + dt3) / 6.0;
        let p = FourMomentum((1.0 + u.squared()).sqrt(), u.0, u.1, u.2) * m;

        ChargedParticle {
            charge: e,
            position: x,
            momentum: p,
        }
    }
}

// ============================================================================
#[cfg(test)]
mod test {

    use crate::charged_particle::ChargedParticle;
    use crate::four_momentum::FourMomentum;
    use crate::three_vector::ThreeVector;

    #[test]
    fn rk4_push_works_with_zero_fields() {
        let particle = ChargedParticle {
            charge: 1.0,
            position: ThreeVector(0.0, 0.0, 0.0),
            momentum: FourMomentum::from_mass_and_velocity(1.0, ThreeVector(0.5, 0.0, 0.0)),
        };

        let field = |_, _| {
            let e = ThreeVector(0.0, 0.0, 0.0);
            let b = ThreeVector(0.0, 0.0, 0.0);
            (e, b)
        };

        let p1 = particle.rk4_push(field, 0.0, 0.1);

        approx::assert_relative_eq!(p1.position.0, 0.05, epsilon = 1e-12);
        approx::assert_relative_eq!(p1.momentum.0, particle.momentum.0, epsilon = 1e-12);
        approx::assert_relative_eq!(p1.momentum.1, particle.momentum.1, epsilon = 1e-12);
    }
}
