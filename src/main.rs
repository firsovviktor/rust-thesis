use std::sync::Arc;

use tokio::runtime::Builder;

const NUM_THREADS: usize = 2;

#[derive(Default)]
struct Solution {
    black: Vec<f64>,
    red: Vec<f64>,
}

impl Solution {
    fn new_black(idx: usize, total_size: usize, prev_black: f64, prev_red: f64) -> f64 {
        prev_black + prev_red
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

        black_center + prev_red
    }

    async fn calculate_black(
        left: usize,
        right: usize,
        black: Arc<Vec<f64>>,
        red: Arc<Vec<f64>>,
        out_black: Arc<Vec<f64>>,
    ) {
        for i in left..right {
            unsafe {
                // Get a mutable pointer to the vector's data
                let ptr = out_black.as_ptr() as *mut f64;
                // Increment the value at the calculated index
                *ptr.add(i) = Self::new_black(i, black.len(), black[i], red[i])
            }
        }
    }

    async fn calculate_red(
        left: usize,
        right: usize,
        black: Arc<Vec<f64>>,
        red: Arc<Vec<f64>>,
        out_red: Arc<Vec<f64>>,
    ) {
        for i in left..right {
            unsafe {
                // Get a mutable pointer to the vector's data
                let ptr = out_red.as_ptr() as *mut f64;
                // Increment the value at the calculated index
                *ptr.add(i) = Self::new_red(
                    i,
                    red.len(),
                    if i > 0 {
                        black.get(i - 1).copied()
                    } else {
                        None
                    },
                    black[i],
                    black.get(i + 1).copied(),
                    red[i],
                );
            }
        }
    }

    async fn step(
        &self,
        black: Arc<Vec<f64>>,
        red: Arc<Vec<f64>>,
    ) -> (Arc<Vec<f64>>, Arc<Vec<f64>>) {
        let mut black_vec = Vec::new();
        black_vec.resize(black.len(), 0.);

        let mut new_black = Arc::new(black_vec);

        let mut tasks = Vec::new();
        for i in 0..NUM_THREADS {
            let left = i * black.len() / NUM_THREADS;
            let right = (i + 1) * black.len() / NUM_THREADS;

            tasks.push(tokio::spawn(Self::calculate_black(
                left,
                right,
                black.clone(),
                red.clone(),
                new_black.clone(),
            )));
        }

        for task in tasks {
            task.await.unwrap();
        }

        let mut red_vec = Vec::new();
        red_vec.resize(red.len(), 0.);

        let mut new_red = Arc::new(red_vec);

        let mut tasks = Vec::new();

        for i in 0..NUM_THREADS {
            let left = i * red.len() / NUM_THREADS;
            let right = (i + 1) * red.len() / NUM_THREADS;

            tasks.push(tokio::spawn(Self::calculate_red(
                left,
                right,
                black.clone(),
                red.clone(),
                new_red.clone(),
            )));
        }

        for task in tasks {
            task.await.unwrap();
        }

        (new_black, new_red)
    }

    async fn calculate(&self) {
        let mut black = Arc::new(vec![1., 1., 1., 1., 1., 1., 1.]);
        let mut red = Arc::new(vec![1., 1., 1., 1., 1., 1., 1.]);

        for _ in 0..10 {
            let (new_black, new_red) = self.step(black, red).await;
            black = new_black;
            red = new_red;
        }

        println!("{:?}", black);
        println!("{:?}", red);
        // println!("{}", black[0]);
    }
}

fn main() {
    let runtime = Builder::new_multi_thread()
        .worker_threads(NUM_THREADS)
        .enable_all()
        .build()
        .unwrap();

    let solution = Solution::default();

    runtime.block_on(solution.calculate());
}
