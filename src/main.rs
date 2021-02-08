#![allow(unused)]

use crate::charged_particle::ChargedParticle;
use crate::four_momentum::FourMomentum;
use crate::three_vector::ThreeVector;
use std::fs::File;
use std::io::prelude::*;
use std::io::Write;

mod charged_particle;
mod four_momentum;
mod three_vector;

fn main() {
    let mut particle = ChargedParticle {
        charge: 1.0,
        position: ThreeVector(0.0, 0.0, 0.0),
        momentum: FourMomentum::from_mass_and_velocity(1.0, ThreeVector(0.0, 0.5, 0.0)),
    };
    let mut time: f64 = 0.0;
    let dt: f64 = 0.1;

    let field = |_, _| {
        let e = ThreeVector(0.0, 0.2, 0.0);
        let b = ThreeVector(0.0, 0.0, 1.0);
        (e, b)
    };

    let save_particle_vec = |solution: &Vec<ChargedParticle>| {
        let mut file = File::create("solution.dat").expect("error creating");
        for element in solution {
            let x = element.position.0;
            let y = element.position.1;
            let z = element.position.2;
            writeln!(&mut file, "{} {} {}", x, y, z);
        }
    };
    let mut solution = vec![];
    //let e = ThreeVector(0.0, 0.2, 0.0);
    //let b = ThreeVector(0.0, 0.0, 1.0);
    //
    while time < 16.0 {
        particle = particle.rk4_push(field, time, dt);
        //      println!(
        //          "{} {} {} {}",
        //          time, particle.position.0, particle.position.1, particle.position.2
        //      );

        solution.push(particle.clone());

        time += dt;
    }
    //println!("{:?}", solution);
    save_particle_vec(&solution);
}
