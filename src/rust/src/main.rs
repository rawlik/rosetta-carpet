use plotters::prelude::*;
use std::ops::Range;
use std::f64::consts::PI;
use std::sync::Arc;
use itertools_num::linspace;
use rustfft::num_complex::Complex;
use rustfft::num_traits::Zero;
use rustfft::num_traits::One;
use rustfft::FFTplanner;
use ndarray::Array;


fn fftfreq(n: usize, d: f64) -> Vec<f64> {
    // ref. https://docs.scipy.org/doc/numpy/reference/generated/numpy.fft.fftfreq.html
    let parity = n % 2;
    let n = n as i32;
    match parity {
        // f = [0, 1, ...,   n/2-1,     -n/2, ..., -1] / (d*n)   if n is even
        0 => (0..n/2).chain(-n/2..0),
        // f = [0, 1, ..., (n-1)/2, -(n-1)/2, ..., -1] / (d*n)   if n is odd
        1 => (0..(n-1)/2+1).chain(-(n-1)/2..0),
        _ => panic!("i32 % 2 was neither 0 nor 1!")
    }
        .map(|f| f as f64 / (d * n as f64))
        .collect()
}


fn main() -> Result<(), Box<dyn std::error::Error>> {
    let period: f64 = 10e-6;
    let lambda = 26.95e-12;
    let dy = 0.1e-6;
    let talbot_distance = 2. * period.powi(2) / lambda;
    let ny = 600;
    let nz = 800;
    // let dz = 0.1e-6;
    let dz = talbot_distance / nz as f64;
    println!("The talbot_distance is {:?}m", talbot_distance);
    let fy = fftfreq(ny, dy);

    // The initial wavefront
    // let mut y = Array::range(0., (ny as f64) * dy, dy)
    let y: Vec<Complex<f64>> = (0..ny)
        .map(|y| (y as f64) * dy)
        .map(|y| 0.5 * (1. + (2. * PI * y / period).cos().signum()))
        .map(|y| Complex::new(y, 0.))
        .collect();
    // println!("{:?}", y);

    let mut planner = FFTplanner::new(false);
    let fft = planner.plan_fft(ny);
    let mut iplanner = FFTplanner::<f64>::new(true);
    let ifft = iplanner.plan_fft(ny);
    let mut input: Vec<Complex<f64>> = vec![Complex::zero(); ny];
    let mut output: Vec<Complex<f64>> = vec![Complex::zero(); ny];

    let mut z: Vec<Vec<Complex<f64>>> = vec![vec![Complex::zero(); ny]; nz];

    z[0].copy_from_slice(&y);

    for iz in 1..nz {
        // fft
        input.copy_from_slice(&z[iz - 1]);
        fft.process(&mut input, &mut z[iz]);

        // processing in the Fourier space
        for (v, f) in z[iz].iter_mut().zip(fy[..ny].iter()) {
            // normalization
            *v *= Complex::new(1.0 / (ny as f64), 0.0);
            // transfer function for the Fresnel free propagation
            let ddz = dz;// * iz as f64;
            *v *= Complex::new(0., 2.0 * PI / lambda * ddz).exp();
            *v *= Complex::new(0., -1. * PI * lambda * ddz * f.powi(2)).exp();
        }

        // inverse fft
        input.copy_from_slice(&z[iz]);
        ifft.process(&mut input, &mut z[iz]);
    }


    // fft.process(&mut y.clone(), &mut output);

    // The fft instance returned by the planner is stored behind an `Arc`, so it's cheap to clone
    // let fft_clone = Arc::clone(&fft);


    let image_path = std::env::current_exe()?
        .parent().expect("Wrong path.")
        .parent().expect("Wrong path.")
        .parent().expect("Wrong path.")
        .parent().expect("Wrong path.")
        .parent().expect("Wrong path.")
        .join("docs")
        .join("carpets")
        .join("rust.png");
    println!("{:?}", image_path);

    let root =
        BitMapBackend::new(&image_path, (800, 600)).into_drawing_area();

    root.fill(&White)?;

    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .x_label_area_size(10)
        .y_label_area_size(10)
        .build_ranged(0..800, 0..600)?;

    chart
        .configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        .draw()?;

    let plotting_area = chart.plotting_area();

    let range = plotting_area.get_pixel_range();

    let (pw, ph) = (range.0.end - range.0.start, range.1.end - range.1.start);
    let (xr, yr) = (chart.x_range(), chart.y_range());

    // for (x, y, c) in mandelbrot_set(xr, yr, (pw as usize, ph as usize), 100) {
    //     if c != 100 {
    //         plotting_area.draw_pixel((x, y), &HSLColor(c as f64 / 100.0, 1.0, 0.5))?;
    //     } else {
    //         plotting_area.draw_pixel((x, y), &Black)?;
    //     }
    // }

    // let singrating(x: Iterator, L: f64) = x.map()

    for xplot in 0..800 {
        for yplot in 0..600 {
            plotting_area.draw_pixel(
                (xplot as i32, yplot as i32),
                &HSLColor(0.0, 0.0, z[xplot][yplot].norm()))?;
        }
    }

    // let x: Vec<f64> = (0..Nx).map(|x| (x as f64) * dx - (Nx as f64) / 2.0)
    //                          .collect();
    // let x: Vec<Complex<f64>> = vec![Complex::one()];
    // println!("{:?}", x);
    // let xx: Vec<f64> = x.iter()
    //                      .map(|x| 0.5 * (1.0 + (2.0 * PI * x / 10e-6).cos().signum()))
    //                      .collect();
    // println!("{:?}", xx);

    // let xxx: Vec<Complex<f64>> = xx;

    // let y: i64 = 5;
    // let _z: f64 = y as f64;
    Ok(())
}

fn mandelbrot_set(
    real: Range<f64>,
    complex: Range<f64>,
    samples: (usize, usize),
    max_iter: usize,
) -> impl Iterator<Item = (f64, f64, usize)> {
    let step = (
        (real.end - real.start) / samples.0 as f64,
        (complex.end - complex.start) / samples.1 as f64,
    );
    return (0..(samples.0 * samples.1)).map(move |k| {
        let c = (
            real.start + step.0 * (k % samples.0) as f64,
            complex.start + step.1 * (k / samples.0) as f64,
        );
        let mut z = (0.0, 0.0);
        let mut cnt = 0;
        while cnt < max_iter && z.0 * z.0 + z.1 * z.1 <= 1e10 {
            z = (z.0 * z.0 - z.1 * z.1 + c.0, 2.0 * z.0 * z.1 + c.1);
            cnt += 1;
        }
        return (c.0, c.1, cnt);
    });
}

// fn rect(x: f64) -> {};
// rect(x) = one(x) .* (abs.(x) .< 0.5)