use std::{
    sync::Arc,
    time::{Duration, Instant},
};

use tokio::runtime::Builder;

const NUM_THREADS: usize = 4;

struct ArgsForBlack {
    pub left: usize,
    pub right: usize,
    pub len: usize,
    pub black: *const f64,
    pub red: *const f64,
    pub out_black: *mut f64,
}
unsafe impl Send for ArgsForBlack {}

struct ArgsForRed {
    pub left: usize,
    pub right: usize,
    pub len: usize,
    pub black: *const f64,
    pub red: *const f64,
    pub out_red: *mut f64,
}

unsafe impl Send for ArgsForRed {}

#[derive(Default)]
struct Solution {
    black: Vec<f64>,
    red: Vec<f64>,
    new_black: Vec<f64>,
    new_red: Vec<f64>,
}

impl Solution {
    fn new_black(idx: usize, total_size: usize, prev_black: f64, prev_red: f64) -> f64 {
        prev_black + prev_red % 2.
    }

    fn new_red(
        idx: usize,
        total_size: usize,
        black_left: Option<f64>,
        black_center: f64,
        black_right: Option<f64>,
        prev_red: f64,
    ) -> f64 {
        // let total = black_left.unwrap_or(0.) + black_center + black_right.unwrap_or(0.);
        // let total = black_left.unwrap() + black_center + black_right.unwrap();

        prev_red + black_center % 2.
    }

    async fn calculate_black(args: ArgsForBlack) {
        for i in args.left..args.right {
            unsafe {
                // Get a mutable pointer to the vector's data
                // Increment the value at the calculated index
                *args.out_black.add(i) =
                    Self::new_black(i, args.len, *args.black.add(i), *args.red.add(i))
            }
        }
    }

    async fn calculate_red(args: ArgsForRed) {
        for i in args.left..args.right {
            unsafe {
                // Get a mutable pointer to the vector's data
                // Increment the value at the calculated index
                *args.out_red.add(i) = Self::new_red(
                    i,
                    args.len,
                    if i > 0 {
                        Some(*args.black.add(i - 1))
                    } else {
                        None
                    },
                    *args.black.add(i),
                    if i < args.len - 1 {
                        Some(*args.black.add(i + 1))
                    } else {
                        None
                    },
                    *args.red.add(i),
                );
            }
        }
    }

    async fn step(&mut self) {
        let mut tasks = Vec::new();
        for i in 0..NUM_THREADS {
            let left = i * self.black.len() / NUM_THREADS;
            let right = ((i + 1) * self.black.len() / NUM_THREADS).min(self.black.len());

            tasks.push(tokio::spawn(Self::calculate_black(ArgsForBlack {
                left,
                right,
                len: self.black.len(),
                black: self.black.as_ptr(),
                red: self.red.as_ptr(),
                out_black: self.new_black.as_ptr() as *mut f64,
            })));
        }

        for task in tasks {
            task.await.unwrap();
        }

        let mut tasks = Vec::new();

        for i in 0..NUM_THREADS {
            let left = i * self.red.len() / NUM_THREADS;
            let right = ((i + 1) * self.red.len() / NUM_THREADS).min(self.red.len());

            tasks.push(tokio::spawn(Self::calculate_red(ArgsForRed {
                left,
                right,
                len: self.red.len(),
                black: self.black.as_ptr(),
                red: self.red.as_ptr(),
                out_red: self.new_red.as_ptr() as *mut f64,
            })));
        }

        for task in tasks {
            task.await.unwrap();
        }

        std::mem::swap(&mut self.black, &mut self.new_black);
        std::mem::swap(&mut self.red, &mut self.new_red);
    }

    async fn calculate(&mut self) {
        let t = Instant::now();

        self.black = vec![1.; 100_000];
        self.red = vec![1.; 100_000];
        self.new_black.resize(self.black.len(), 0.);
        self.new_red.resize(self.red.len(), 0.);

        for _ in 0..100000 {
            self.step().await;
        }

        // println!("{:?}", self.black);
        // println!("{:?}", self.red);
        // println!("{}", black[0]);

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

    let mut solution = Solution::default();

    runtime.block_on(solution.calculate());
}
