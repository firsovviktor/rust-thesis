use std::{
    fmt::format,
    fs::File,
    io::{BufWriter, Read, Write},
    time::{Duration, Instant},
};

use tokio::{
    runtime::Builder,
    sync::mpsc::{Receiver, Sender},
};

use num_complex::Complex;
use statrs::statistics::Statistics;

const GLOBAL_NUM_THREADS: usize = 14;
const EXP_THREADS_SPLIT: usize = 2;
const EXP_ENERGY_THREADS_SPLIT: usize = 7;
const OUTPUT_FILES_DIR: &str = "output_files";

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

    pub energy_pulling_frequency: usize,

    pub base_Nr: usize,
    pub base_Nt: usize,
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

struct ArgsForNonLinearity {
    pub left: usize,
    pub right: usize,
    pub field: *const Complex<f64>,
    pub params: *const Parameters,
}
unsafe impl Send for ArgsForNonLinearity {}

struct ArgsForEnergy {
    pub left: usize,
    pub right: usize,
    pub field: *const Complex<f64>,
    pub impulse: *const Complex<f64>,
    pub next_field: *const Complex<f64>,
    pub params: *const Parameters,
}
unsafe impl Send for ArgsForEnergy {}

struct Physics {
    pub Action: Complex<f64>,
    pub Energy: Vec<Complex<f64>>,
}
unsafe impl Send for Physics {}

struct Solution {
    field: Vec<Complex<f64>>,
    next_field: Vec<Complex<f64>>,

    impulse: Vec<Complex<f64>>,
    next_impulse: Vec<Complex<f64>>,

    physics: Physics,

    params: Parameters,
    sender: Sender<DumpElement>,

    output_file: BufWriter<File>,
}

impl Solution {
    // constructor
    pub fn new(params: Parameters, sender: Sender<DumpElement>) -> Self {
        let next_field = vec![Complex::new(0.0, 0.0); params.Nr];
        let next_impulse = vec![Complex::new(0.0, 0.0); params.Nr];
        let (field, impulse) = Solution::initialize_conditions(&params);
        let physics = Physics {
            Action: Complex::new(0.0, 0.0),
            Energy: vec![],
        };

        let out_filename = format!(
            "{}/Nr_{}_Nt_{}_alpha_{}_Dimensions_{}.txt",
            OUTPUT_FILES_DIR, params.Nr, params.Nt, params.alpha, params.Dimensions
        );
        let mut output_file = BufWriter::new(File::create(out_filename).unwrap());
        writeln!(
            &mut output_file,
            "t r field_re field_im impulse_re impulse_im"
        )
        .unwrap();

        Solution {
            field,
            next_field,
            impulse,
            next_impulse,
            physics,
            params,
            sender,
            output_file,
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

    // function calculates non-linearity at given spatial index
    fn new_NonLinearity(index: usize, field: Complex<f64>, params: &Parameters) -> Complex<f64> {
        if index == 0 || index == params.Nr - 1 {
            return params.Lambda * field.powc(Complex::new(4.0, 0.0)) * params.dr / 2.0
                * (params.R0 + (index as f64) * params.dr).powf(params.Dimensions - 1.0);
        } else {
            return params.Lambda
                * field.powc(Complex::new(4.0, 0.0))
                * params.dr
                * (params.R0 + (index as f64) * params.dr).powf(params.Dimensions - 1.0);
        }
    }

    // function calculates energy density at given spatial index
    fn new_Energy(
        index: usize,
        field: Complex<f64>,
        field_right: Option<Complex<f64>>,
        next_field: Complex<f64>,
        next_field_right: Option<Complex<f64>>,
        impulse: Complex<f64>,
        params: &Parameters,
    ) -> Complex<f64> {
        let mut out_Energy: Complex<f64>;
        out_Energy = Complex::new(0.0, 0.0);

        let field_right = field_right.unwrap_or(Complex::new(0.0, 0.0));
        let next_field_right = next_field_right.unwrap_or(Complex::new(0.0, 0.0));

        if index == 0 || index == params.Nr - 1 {
            out_Energy += (0.5 * Complex::<f64>::powc(impulse, Complex::new(2.0, 0.0))
                + 0.5
                    * f64::powf(params.Mass, 2.0)
                    * Complex::<f64>::powc(field + next_field, Complex::new(2.0, 0.0))
                    / 4.0
                + 0.25
                    * params.Lambda
                    * Complex::<f64>::powc(field + next_field, Complex::new(4.0, 0.0))
                    / 16.0)
                * params.dr
                / 2.0
                * f64::powf(
                    params.R0 + (index as f64) * params.dr,
                    params.Dimensions - 1.0,
                );
        } else {
            out_Energy += (0.5 * Complex::<f64>::powc(impulse, Complex::new(2.0, 0.0))
                + 0.5
                    * f64::powf(params.Mass, 2.0)
                    * Complex::<f64>::powc(field + next_field, Complex::new(2.0, 0.0))
                    / 4.0
                + 0.25
                    * params.Lambda
                    * Complex::<f64>::powc(field + next_field, Complex::new(4.0, 0.0))
                    / 16.0)
                * params.dr
                * f64::powf(
                    params.R0 + (index as f64) * params.dr,
                    params.Dimensions - 1.0,
                );
        }

        if index != params.Nr - 1 {
            out_Energy += (Complex::<f64>::powc(field_right - field, Complex::new(2.0, 0.0))
                / 4.0
                / params.dr
                + Complex::<f64>::powc(next_field_right - next_field, Complex::new(2.0, 0.0))
                    / 4.0
                    / params.dr)
                * f64::powf(
                    params.R0 + (index as f64) * params.dr + params.dr / 2.0,
                    params.Dimensions - 1.0,
                );
        }

        out_Energy
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

    // function calculates non-linearity for the current field row
    async fn calculate_NonLinearity(args: ArgsForNonLinearity) -> Complex<f64> {
        let mut res: Complex<f64>;
        res = Complex::new(0.0, 0.0);

        for i in args.left..args.right {
            unsafe {
                let params = args.params.as_ref().unwrap();
                res += Self::new_NonLinearity(i, *args.field.add(i), params);
            }
        }
        res
    }

    // function calculates energy for the current field row
    async fn calculate_Energy(args: ArgsForEnergy) -> Complex<f64> {
        let mut res: Complex<f64>;
        res = Complex::new(0.0, 0.0);

        for i in args.left..args.right {
            unsafe {
                let params = args.params.as_ref().unwrap();
                res += Self::new_Energy(
                    i,
                    *args.field.add(i),
                    if i < params.Nr - 1 {
                        Some(*args.field.add(i + 1))
                    } else {
                        None
                    },
                    *args.next_field.add(i),
                    if i < params.Nr - 1 {
                        Some(*args.next_field.add(i + 1))
                    } else {
                        None
                    },
                    *args.impulse.add(i),
                    params,
                );
            }
        }
        res
    }

    // function does one time step, returns vectors of new impulse and field rows
    async fn step(&mut self, time_index: usize) {
        let mut tasks = Vec::new();
        for i in 0..EXP_THREADS_SPLIT {
            let left = i * self.params.Nr / EXP_THREADS_SPLIT;
            let right = ((i + 1) * self.params.Nr / EXP_THREADS_SPLIT).min(self.params.Nr);

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
        for i in 0..EXP_THREADS_SPLIT {
            let left = i * self.params.Nr / EXP_THREADS_SPLIT;
            let right = ((i + 1) * self.params.Nr / EXP_THREADS_SPLIT).min(self.params.Nr);

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

        let mut NonLinearity: Complex<f64>;
        NonLinearity = Complex::new(0.0, 0.0);

        let mut tasks = Vec::new();
        for i in 0..EXP_THREADS_SPLIT {
            let left = i * self.params.Nr / EXP_THREADS_SPLIT;
            let right = ((i + 1) * self.params.Nr / EXP_THREADS_SPLIT).min(self.params.Nr);

            tasks.push(tokio::spawn(Self::calculate_NonLinearity(
                ArgsForNonLinearity {
                    left,
                    right,
                    field: self.field.as_ptr(),
                    params: &self.params as *const Parameters,
                },
            )));
        }

        for task in tasks {
            let result = task.await.unwrap();
            NonLinearity += result;
        }

        //Extra division by 2 is caused by Im S ~ 1/2 Nl ~ 1/2 \int \lambda field^4
        if time_index == 0 || time_index == self.params.Nt - 1 {
            self.physics.Action += NonLinearity / 4.0 * self.params.dt;
        } else {
            self.physics.Action += NonLinearity / 2.0 * self.params.dt;
        }
        if time_index % self.params.energy_pulling_frequency == 0 {
            let mut current_energy: Complex<f64>;
            current_energy = Complex::new(0.0, 0.0);

            let mut tasks = Vec::new();
            for i in 0..EXP_ENERGY_THREADS_SPLIT {
                let left = i * self.params.Nr / EXP_ENERGY_THREADS_SPLIT;
                let right =
                    ((i + 1) * self.params.Nr / EXP_ENERGY_THREADS_SPLIT).min(self.params.Nr);

                tasks.push(tokio::spawn(Self::calculate_Energy(ArgsForEnergy {
                    left,
                    right,
                    field: self.field.as_ptr(),
                    impulse: self.impulse.as_ptr(),
                    next_field: self.next_field.as_ptr(),
                    params: &self.params as *const Parameters,
                })));
            }

            for task in tasks {
                let result = task.await.unwrap();
                current_energy += result;
            }
            self.physics.Energy.push(current_energy);
        }

        let scale: usize = (self.params.Nt - 1) / (self.params.base_Nt - 1);
        //if self.params.Nr == 33 {
        //    println!("sc {}, Nr {}, Nt {}", scale, self.params.Nr, self.params.Nt);
        //}
        if time_index % (scale) == 0 {
            for i in 0..self.params.base_Nr {
                writeln!(
                    &mut self.output_file,
                    "{} {} {} {} {} {}",
                    self.params.T0 + (time_index as f64) * self.params.dt,
                    self.params.R0 + ((i * scale) as f64) * self.params.dr,
                    self.field[i * scale].re,
                    self.field[i * scale].im,
                    self.impulse[i * scale].re,
                    self.impulse[i * scale].im
                )
                .unwrap();
            }
        }

        std::mem::swap(&mut self.field, &mut self.next_field);
        std::mem::swap(&mut self.impulse, &mut self.next_impulse);
    }

    // function repeatedly calucaltes new rows and stores them in place of the old ones
    // also gives opprtunity to caltucalte signatures and store values
    async fn calculate(mut self) {
        let t = Instant::now();

        for i in 0..self.params.Nt {
            self.step(i).await;
        }

        let mean_energy = Complex::new(
            self.physics.Energy.iter().map(|c| c.re).mean(),
            self.physics.Energy.iter().map(|c| c.im).mean(),
        );

        self.sender
            .send(DumpElement {
                task_id: format!(
                    "Nr_{}_Nt_{}_alpha_{}_Dimensions_{}",
                    self.params.Nr, self.params.Nt, self.params.alpha, self.params.Dimensions
                ),
                Action: self.physics.Action,
                Energy: mean_energy,
                Energy_std: (self
                    .physics
                    .Energy
                    .iter()
                    .map(|c| (c.re - mean_energy.re).powf(2.0) + (c.im - mean_energy.im).powf(2.0))
                    .sum::<f64>()
                    / (self.physics.Energy.len() as f64 - 1.0))
                    .powf(0.5),
                elapsed: Some(t.elapsed()),
                terminate: false,
            })
            .await
            .unwrap();

        //println!("{:?}", self.field);
        // println!("Thread elapsed: {}s", t.elapsed().as_secs_f64());
    }
}
unsafe impl Send for Solution {}

fn build_solution(mut params: Parameters, sender: Sender<DumpElement>) -> Solution {
    params.dr = (params.Lr - params.R0) / (params.Nr as f64 - 1.0);
    params.dt = (params.Lt - params.T0) / (params.Nt as f64 - 1.0);

    Solution::new(params, sender)
}

struct DumpElement {
    task_id: String,
    Action: Complex<f64>,
    Energy: Complex<f64>,
    Energy_std: f64,
    elapsed: Option<Duration>,
    terminate: bool,
}
impl DumpElement {
    fn terminate() -> Self {
        Self {
            task_id: "0".to_string(),
            Action: Complex::i(),
            Energy: Complex::i(),
            Energy_std: 0.0,
            elapsed: None,
            terminate: true,
        }
    }
}

async fn info_writer(mut receiver: Receiver<DumpElement>, total_solutions_count: usize) {
    let mut output_file = BufWriter::new(File::create("output.txt").unwrap());

    let mut completed_solutions = 0;

    while let Some(info) = receiver.recv().await {
        if info.terminate {
            return;
        }

        writeln!(
            &mut output_file,
            "task_id: {}, Action_re: {}, Action_im: {}, Energy_re: {}, Energy_im: {}, Energy_std: {}",
            info.task_id, info.Action.re, info.Action.im, info.Energy.re, info.Energy.im, info.Energy_std
        )
        .unwrap();

        completed_solutions += 1;
        println!(
            "Completed: {}/{} ({}%). Completed task: {}. Elapsed: {}s",
            completed_solutions,
            total_solutions_count,
            (completed_solutions * 100) as f64 / total_solutions_count as f64,
            info.task_id,
            info.elapsed.unwrap().as_secs_f64()
        );
    }
}

fn main() {
    let runtime = Builder::new_multi_thread()
        .worker_threads(GLOBAL_NUM_THREADS)
        .build()
        .unwrap();

    std::fs::create_dir_all(OUTPUT_FILES_DIR).unwrap();

    let mut params = vec![];
    let mut sample_param = Parameters {
        Nr: 65,
        Nt: 257,

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

        energy_pulling_frequency: 1,

        base_Nr: 33,
        base_Nt: 129,
    };

    let (sender, receiver) = tokio::sync::mpsc::channel(10000);

    let vec_Nr = [
        (sample_param.base_Nr as i32 - 1) * 2i32.pow(0) + 1,
        (sample_param.base_Nr as i32 - 1) * 2i32.pow(1) + 1,
        (sample_param.base_Nr as i32 - 1) * 2i32.pow(2) + 1,
        (sample_param.base_Nr as i32 - 1) * 2i32.pow(3) + 1,
        (sample_param.base_Nr as i32 - 1) * 2i32.pow(4) + 1,
        (sample_param.base_Nr as i32 - 1) * 2i32.pow(5) + 1,
        (sample_param.base_Nr as i32 - 1) * 2i32.pow(6) + 1,
        (sample_param.base_Nr as i32 - 1) * 2i32.pow(7) + 1,
        (sample_param.base_Nr as i32 - 1) * 2i32.pow(8) + 1,
        (sample_param.base_Nr as i32 - 1) * 2i32.pow(9) + 1,
        (sample_param.base_Nr as i32 - 1) * 2i32.pow(10) + 1,
    ];
    let vec_alpha = [
        -2.5, -2.0, -1.5, -1.2, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0,
        2.5, 2.0, 1.5, 1.2, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1,
    ];
    let vec_Dimensions = [
        0.0, 0.5, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0,
        4.5, 5.0, 6.0, 8.0,
    ];

    //expecting 7018+1 files 750 kB each
    //total size ~5.2 GB
    //total expected worktime on a 8C/16T CPU ~16.7 hours

    for Dimensions in vec_Dimensions {
        for alpha in vec_alpha {
            for Nr in vec_Nr {
                sample_param.Nr = Nr as usize;
                sample_param.Nt = (4 * (Nr - 1) + 1) as usize;
                sample_param.alpha = alpha;
                sample_param.Dimensions = Dimensions;
                if Nr > 1000 {
                    sample_param.energy_pulling_frequency = 10;
                }
                params.push(sample_param);
            }
        }
    }

    println!("Total tasks: {}", params.len());

    let mut solutions = params
        .into_iter()
        .map(|p| build_solution(p, sender.clone()))
        .collect::<Vec<_>>();

    let solutions_len = solutions.len();

    let t = Instant::now();

    runtime.block_on(async move {
        futures::join!(
            async {
                futures::future::join_all(
                    solutions.into_iter().map(|s| tokio::spawn(s.calculate())),
                )
                .await;

                sender.send(DumpElement::terminate()).await.unwrap();
            },
            info_writer(receiver, solutions_len),
        );
    });

    println!("Total elapsed: {}", t.elapsed().as_secs_f64());
}
