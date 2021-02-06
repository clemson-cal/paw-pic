#![allow(unused)]

use crate::charged_particle::ChargedParticle;

use crate::three_vector::ThreeVector;

use crate::four_momentum::FourMomentum;

mod charged_particle;
mod four_momentum;
mod three_vector;

fn main() {
    let particle = ChargedParticle {
        charge: 1.0,
        position: ThreeVector(0.0, 0.0, 0.0),
        momentum: FourMomentum::from_mass_and_velocity(1.0, ThreeVector(0.5, 0.0, 0.0)),
    };
    //let a = three_vector::ThreeVector(1.0, 2.0, 3.0);
    //let b = three_vector::ThreeVector(1.0, 2.0, 3.0);
    let mut time: f64 = 0.0;
    let dt: f64 = 0.1;

    //let field = |_, _| {
    //    let e = ThreeVector(0.0, 0.0, 0.0);
    //    let b = ThreeVector(0.0, 0.0, 1.0);
    //    (e, b)
    //};
    let e = ThreeVector(0.0, 0.0, 0.0);
    let b = ThreeVector(0.0, 0.0, 1.0);

    while time < 16.0 {
        let p1 = particle.boris_push(e, b, dt);
        println!(
            "{} {} {} {}",
            time, p1.position.0, p1.position.1, p1.position.2
        );
        time += dt;
    }

    //println!("dot(a, b) = {}", a.dot(&b));
}
