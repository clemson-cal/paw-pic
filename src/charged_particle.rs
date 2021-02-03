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

        let u_nmh     = self.momentum.gamma_beta_vector();
        let u_minus   = u_nmh + electric * (0.5 * h);
        let gamma_n   = (1.0 + u_minus.squared()).sqrt();
        let t         = magnetic * (0.5 * h / gamma_n);
        let s         = t / (1.0 + t.squared()) * 2.0;
        let u_prime   = u_minus + u_minus.cross(&t);
        let u_plus    = u_minus + u_prime.cross(&s);
        let u_nph     = u_plus + electric * (0.5 * h);
        let gamma_nph = (1.0 + u_nph.squared()).sqrt();
        let new_pos   = self.position + u_nph * dt / gamma_nph;
        let new_mom   = FourMomentum(gamma_nph, u_nph.0, u_nph.1, u_nph.2) * m;

        ChargedParticle {
            charge: e,
            position: new_pos,
            momentum: new_mom,
        }
    }
}
