use std::{
    sync::Arc,
    time::{Duration, Instant},
};

use tokio::runtime::Builder;

use num_complex::Complex;

const NUM_THREADS: usize = 2;

#[derive(Copy, Clone)]
struct Parameters {
    pub Nr: usize,
    pub Nt: usize,

    pub R0: f64,
    pub Lr: f64,
    pub T0: f64,
    pub Lt: f64,

    pub Mass: f64,
    pub Lambda: f64,
    pub Dimensions: f64,

    pub alpha: f64,

    pub r1: f64,
    pub a1: f64,
    pub r2: f64,
    pub a2: f64,
    pub k: f64,

    pub ampl: f64,

    pub dr: f64,
    pub dt: f64,
}

struct ArgsForField {
    pub left: usize,
    pub right: usize,
    pub field: *const Complex<f64>,
    pub impulse: *const Complex<f64>,
    pub out_field: *mut Complex<f64>,
    pub params: *const Parameters,
}
unsafe impl Send for ArgsForField {}

struct ArgsForImpulse {
    pub left: usize,
    pub right: usize,
    pub field: *const Complex<f64>,
    pub impulse: *const Complex<f64>,
    pub out_impulse: *mut Complex<f64>,
    pub params: *const Parameters,
}
unsafe impl Send for ArgsForImpulse {}

struct Solution {
    field: Vec<Complex<f64>>,
    next_field: Vec<Complex<f64>>,

    impulse: Vec<Complex<f64>>,
    next_impulse: Vec<Complex<f64>>,

    params: Parameters,
}

impl Solution {
    // constructor
    pub fn new(params: Parameters) -> Self {
        let next_field = vec![Complex::new(0.0, 0.0); params.Nr];
        let next_impulse = vec![Complex::new(0.0, 0.0); params.Nr];
        let (field, impulse) = Solution::initialize_conditions(&params);
        Solution {
            field,
            next_field,
            impulse,
            next_impulse,
            params,
        }
    }

    // function initializes field and impulse at the strating time
    fn initialize_conditions(params: &Parameters) -> (Vec<Complex<f64>>, Vec<Complex<f64>>) {
        let mut field = vec![Complex::new(0.0, 0.0); params.Nr];
        let mut impulse = vec![Complex::new(0.0, 0.0); params.Nr];

        let mut r: f64;
        for i in 0..params.Nr {
            r = params.R0 + params.dr * i as f64 + params.T0;
            field[i] = Complex::new(0.0, 1.0) * params.alpha * params.ampl * f64::cos(params.k * r)
                + params.ampl * f64::sin(params.k * r);
            field[i] /= 1.0 + f64::exp(-params.a1 * (r - params.r1));
            field[i] /= 1.0 + f64::exp(-params.a2 * (params.r2 - r));
            if i > 0 {
                field[i] /= f64::powf(params.dr * i as f64, (params.Dimensions - 1.0) / 2.0);
            }
        }

        let mut pi1: Complex<f64>;
        let mut pi2: Complex<f64>;
        let mut pi3: Complex<f64>;

        for i in 0..params.Nr {
            r = params.R0 + params.dr * i as f64 + params.T0 + params.dt / 2.0;

            pi1 = (Complex::new(0.0, 1.0) * params.alpha * params.ampl * f64::cos(params.k * r)
                + params.ampl * f64::sin(params.k * r))
                * params.a1;
            pi1 *= f64::exp(-params.a1 * (r - params.r1));
            pi1 /= f64::powf(1.0 + f64::exp(-params.a1 * (r - params.r1)), 2.0);
            pi1 /= 1.0 + f64::exp(-params.a2 * (params.r2 - r));

            pi2 = -(Complex::new(0.0, 1.0) * params.alpha * params.ampl * f64::cos(params.k * r)
                + params.ampl * f64::sin(params.k * r))
                * params.a2;
            pi2 *= f64::exp(-params.a2 * (params.r2 - r));
            pi2 /= 1.0 + f64::exp(-params.a1 * (r - params.r1));
            pi2 /= f64::powf(1.0 + f64::exp(-params.a2 * (params.r2 - r)), 2.0);

            pi3 = -Complex::new(0.0, 1.0)
                * params.alpha
                * params.ampl
                * params.k
                * f64::sin(params.k * r)
                + params.ampl * params.k * f64::cos(params.k * r);
            pi3 /= 1.0 + f64::exp(-params.a1 * (r - params.r1));
            pi3 /= 1.0 + f64::exp(-params.a2 * (params.r2 - r));

            impulse[i] = pi1 + pi2 + pi3;
            if i > 0 {
                impulse[i] /= f64::powf(params.dr * i as f64, (params.Dimensions - 1.0) / 2.0);
            }
        }
        (field, impulse)
    }

    // function calculates next field value
    fn new_field(field: Complex<f64>, impulse: Complex<f64>, params: &Parameters) -> Complex<f64> {
        field + impulse * params.dt
    }

    // function calculates next impulse value
    fn new_impulse(
        index: usize,
        field_left: Option<Complex<f64>>,
        field_center: Complex<f64>,
        field_right: Option<Complex<f64>>,
        impulse: Complex<f64>,
        params: &Parameters,
    ) -> Complex<f64> {
        let field_left = field_left.unwrap_or(Complex::new(0.0, 0.0));
        let field_right = field_right.unwrap_or(Complex::new(0.0, 0.0));

        let mut out_impulse: Complex<f64>;
        let field_pp: Complex<f64>;
        let field_p: Complex<f64>;

        if index == 0 {
            field_pp = 2.0 * (field_right - field_center) / f64::powf(params.dr, 2.0);
            out_impulse = impulse;
            out_impulse += (params.Dimensions * field_pp
                - f64::powf(params.Mass, 2.0) * field_center)
                * params.dt;
            out_impulse -= params.Lambda * field_center.powc(Complex::new(3.0, 0.0)) * params.dt;
        } else if index == params.Nr - 1 {
            field_pp = 2.0 * (field_left - field_center) / f64::powf(params.dr, 2.0);
            out_impulse = impulse;
            out_impulse += (field_pp - f64::powf(params.Mass, 2.0) * field_center) * params.dt;
            out_impulse -= params.Lambda * field_center.powc(Complex::new(3.0, 0.0)) * params.dt;
        } else {
            field_pp = (field_right - 2.0 * field_center + field_left) / f64::powf(params.dr, 2.0);
            field_p = (field_right - field_left) / 2.0 / params.dr;
            out_impulse = impulse;
            out_impulse += (field_pp
                + (params.Dimensions - 1.0) * field_p / (params.R0 + (index as f64) * params.dr))
                * params.dt;
            out_impulse -= f64::powf(params.Mass, 2.0) * field_center * params.dt;
            out_impulse -= params.Lambda * field_center.powc(Complex::new(3.0, 0.0)) * params.dt;
        }

        out_impulse
    }

    // function calculates next field row
    async fn calculate_field(args: ArgsForField) {
        for i in args.left..args.right {
            unsafe {
                // Get a mutable pointer to the vector's data
                // Increment the value at the calculated index
                *args.out_field.add(i) = Self::new_field(
                    *args.field.add(i),
                    *args.impulse.add(i),
                    args.params.as_ref().unwrap(),
                )
            }
        }
    }

    // function calculates next impulse row
    async fn calculate_impulse(args: ArgsForImpulse) {
        for i in args.left..args.right {
            unsafe {
                let params = args.params.as_ref().unwrap();
                // Get a mutable pointer to the vector's data
                // Increment the value at the calculated index
                *args.out_impulse.add(i) = Self::new_impulse(
                    i,
                    if i > 0 {
                        Some(*args.field.add(i - 1))
                    } else {
                        None
                    },
                    *args.field.add(i),
                    if i < params.Nr - 1 {
                        Some(*args.field.add(i + 1))
                    } else {
                        None
                    },
                    *args.impulse.add(i),
                    params,
                );
            }
        }
    }

    // function does one time step, returns vectors of new impulse and field rows
    async fn step(&mut self) {
        let mut tasks = Vec::new();
        for i in 0..NUM_THREADS {
            let left = i * self.params.Nr / NUM_THREADS;
            let right = ((i + 1) * self.params.Nr / NUM_THREADS).min(self.params.Nr);

            tasks.push(tokio::spawn(Self::calculate_field(ArgsForField {
                left,
                right,
                field: self.field.as_ptr(),
                impulse: self.impulse.as_ptr(),
                out_field: self.next_field.as_ptr() as *mut Complex<f64>,
                params: &self.params as *const Parameters,
            })));
        }

        for task in tasks {
            task.await.unwrap();
        }

        let mut tasks = Vec::new();

        for i in 0..NUM_THREADS {
            let left = i * self.params.Nr / NUM_THREADS;
            let right = ((i + 1) * self.params.Nr / NUM_THREADS).min(self.params.Nr);

            tasks.push(tokio::spawn(Self::calculate_impulse(ArgsForImpulse {
                left,
                right,
                field: self.next_field.as_ptr(),
                impulse: self.impulse.as_ptr(),
                out_impulse: self.next_impulse.as_ptr() as *mut Complex<f64>,
                params: &self.params as *const Parameters,
            })));
        }

        for task in tasks {
            task.await.unwrap();
        }

        std::mem::swap(&mut self.field, &mut self.next_field);
        std::mem::swap(&mut self.impulse, &mut self.next_impulse);
    }

    // function repeatedly calucaltes new rows and stores them in place of the old ones
    // also gives opprtunity to caltucalte signatures and store values
    async fn calculate(&mut self) {
        let t = Instant::now();

        for i in 0..self.params.Nt {
            println!("{}", i);
            println!("{:?}", self.field);
            self.step().await;
        }

        //println!("{:?}", self.field);
        println!("Elapsed: {}s", t.elapsed().as_secs_f64());
    }
}
unsafe impl Send for Solution {}

fn main() {
    let runtime = Builder::new_multi_thread()
        .worker_threads(NUM_THREADS)
        .thread_keep_alive(Duration::from_secs(5))
        .build()
        .unwrap();

    let mut params = Parameters {
        Nr: 17,
        Nt: 65,

        R0: 0.0, //change boundary conditions for non-zero R0
        Lr: 4.0,
        T0: 0.0,
        Lt: 8.0,

        Mass: 1.0,
        Lambda: 1.0,
        Dimensions: 1.0,

        alpha: 0.1,

        r1: 1.5,
        a1: 10.0,
        r2: 2.5,
        a2: 10.0,
        k: 2.0 * std::f64::consts::PI / 1.0,

        ampl: 1.0,

        dr: 0.0,
        dt: 0.0,
    };
    params.dr = (params.Lr - params.R0) / (params.Nr as f64 - 1.0);
    params.dt = (params.Lt - params.T0) / (params.Nt as f64 - 1.0);

    let mut solution = Solution::new(params);

    runtime.block_on(solution.calculate());
}
