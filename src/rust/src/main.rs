use plotters::prelude::*;
use std::convert::TryFrom;
use std::f64::consts::PI;
use rustfft::num_complex::Complex;
use rustfft::num_traits::Zero;
use rustfft::FFTplanner;


fn main() -> Result<(), Box<dyn std::error::Error>> {
    let period: f64 = 10e-6;
    let lambda = 26.95e-12;

    let ny: usize = 600;
    let nz: usize = 800;

    let talbot_distance = 2. * period.powi(2) / lambda;
    let dz = talbot_distance / nz as f64;
    let dy = period * 5. / ny as f64;
    let fy = fftfreq(ny, dy);

    // The initial wavefront - rectangular grating
    let g0: Vec<Complex<f64>> = (0..ny)
        .map(|y| (y as f64) * dy)
        .map(|y| 0.5 * (1. + (2. * PI * y / period).cos().signum()))
        .map(|y| Complex::new(y, 0.))
        .collect();

    // prepare the FFT package, also for the inverse transform
    let mut planner = FFTplanner::new(false);
    let fft = planner.plan_fft(ny);
    let mut iplanner = FFTplanner::<f64>::new(true);
    let ifft = iplanner.plan_fft(ny);
    let mut input: Vec<Complex<f64>> = vec![Complex::zero(); ny];

    // where the carpet is stored
    let mut z: Vec<Vec<Complex<f64>>> = vec![vec![Complex::zero(); ny]; nz];
    z[0].copy_from_slice(&g0);

    for iz in 1..nz {
        // fft
        input.copy_from_slice(&z[0]);
        fft.process(&mut input, &mut z[iz]);

        // processing in the Fourier space
        for (v, f) in z[iz].iter_mut().zip(fy[..ny].iter()) {
            // normalization
            *v *= Complex::new(1.0 / (ny as f64), 0.0);
            // transfer function for the Fresnel free propagation
            let ddz = dz * iz as f64;
            *v *= Complex::new(0., 2.0 * PI / lambda * ddz).exp();
            *v *= Complex::new(0., -1. * PI * lambda * ddz * f.powi(2)).exp();
        }

        // inverse fft
        input.copy_from_slice(&z[iz]);
        ifft.process(&mut input, &mut z[iz]);
    }

    // plotting
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


    let nz = u32::try_from(nz)?;
    let ny = u32::try_from(ny)?;

    let root =
        BitMapBackend::new(&image_path, (nz, ny)).into_drawing_area();

    root.fill(&White)?;

    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .x_label_area_size(10)
        .y_label_area_size(10)
        .build_ranged(0..nz, 0..ny)?;

    chart
        .configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        .draw()?;

    let plotting_area = chart.plotting_area();

    for xplot in 0..nz {
        for yplot in 0..ny {
            plotting_area.draw_pixel(
                (xplot, yplot),
                &HSLColor(0.0, 0.0, z[xplot as usize][yplot as usize].norm()))?;
        }
    }

    Ok(())
}


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
