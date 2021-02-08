use std::ops::{Add, Div, Mul, Sub};

/**
 * Describes a three-dimensional vector
 */
#[derive(Debug, Clone, Copy)]
pub struct ThreeVector(pub f64, pub f64, pub f64);

// ============================================================================
impl ThreeVector {
    pub fn cross(&self, b: &ThreeVector) -> Self {
        ThreeVector(
            self.1 * b.2 - self.2 * b.1,
            self.2 * b.0 - self.0 * b.2,
            self.0 * b.1 - self.1 * b.0,
        )
    }

    pub fn dot(&self, b: &ThreeVector) -> f64 {
        self.0 * b.0 + self.1 * b.1 + self.2 * b.2
    }

    pub fn squared(&self) -> f64 {
        self.dot(self)
    }

    pub fn norm(&self) -> f64 {
        self.squared().sqrt()
    }

    pub fn cosine(&self, b: &ThreeVector) -> f64 {
        self.dot(b) / (self.norm() * b.norm())
    }

    pub fn sine(&self, b: &ThreeVector) -> f64 {
        self.cross(b).norm() / (self.norm() * b.norm())
    }
}

// ============================================================================
impl Add<ThreeVector> for ThreeVector {
    type Output = ThreeVector;
    fn add(self, b: ThreeVector) -> ThreeVector {
        ThreeVector(self.0 + b.0, self.1 + b.1, self.2 + b.2)
    }
}

impl Sub<ThreeVector> for ThreeVector {
    type Output = ThreeVector;
    fn sub(self, b: ThreeVector) -> ThreeVector {
        ThreeVector(self.0 - b.0, self.1 - b.1, self.2 - b.2)
    }
}

impl Mul<f64> for ThreeVector {
    type Output = ThreeVector;
    fn mul(self, b: f64) -> ThreeVector {
        ThreeVector(self.0 * b, self.1 * b, self.2 * b)
    }
}

impl Div<f64> for ThreeVector {
    type Output = ThreeVector;
    fn div(self, b: f64) -> ThreeVector {
        ThreeVector(self.0 / b, self.1 / b, self.2 / b)
    }
}
