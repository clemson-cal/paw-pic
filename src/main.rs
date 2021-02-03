#![allow(unused)]


mod charged_particle;
mod four_momentum;
mod three_vector;


fn main() {

    let a = three_vector::ThreeVector(1.0, 2.0, 3.0);
    let b = three_vector::ThreeVector(1.0, 2.0, 3.0);

    println!("dot(a, b) = {}", a.dot(&b));
}
