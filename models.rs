use ndarray::prelude::*;
use plotters::prelude::*;
use std::f64::consts::PI;

// Define wavelength range (in nm) from 500 to 800, resolution 0.5 nm
let wavelengths: Array1<f64> = Array::range(200.0, 800.5, 0.5);

// Cauchy's equation coefficients (defined for any material)
let a = 1.458;
let b = 0.00354;

// Quantify refractive index (n)
let n_thin_film: Array1<f64> = &wavelengths.mapv(|w| a + b / (w * w));

// Refractive index for base material (constant for simplicity)
let n_si = 3.5;

// Air
let n_air = 1.0;

// Reflectivity using the transfer matrix method to account for interference and multiple reflections
fn transfer_matrix_reflectivity(
    n_air: f64,
    n_thin_film: &Array1<f64>,
    n_silicon: f64,
    thickness: f64,
    wavelengths: &Array1<f64>,
) -> Array1<f64> {
    // Å to nm
    let thickness_nm = thickness * 1e-1;

    // Phase change on reflection
    let delta: Array1<f64> = (2.0 * PI / wavelengths) * n_thin_film * thickness_nm;

    // Reflectivity using transfer matrix method
    let r01: Array1<f64> = (&n_air - n_thin_film) / (&n_air + n_thin_film);
    let r12: Array1<f64> = (n_thin_film - n_silicon) / (n_thin_film + n_silicon);

    let r: Array1<f64> = (&r01 + &r12 * &delta.mapv(|d| (-2.0 * d * 1.0_f64).exp())) /
        (&Array::ones(wavelengths.len()) + &r01 * &r12 * &delta.mapv(|d| (-2.0 * d * 1.0_f64).exp()));

    r.mapv(|r| r.abs().powi(2))
}

// Initialize arrays to store reflectivity spectra for different thicknesses
let thicknesses: Array1<f64> = Array::range(0.0, 6001.0, 1.0); // Thickness range from 0 Å to 6000 Å in steps of 1 Å
let mut reflectivity_spectra: Array2<f64> = Array2::zeros((thicknesses.len(), wavelengths.len()));

// Calculate reflectivity spectra for each thickness
for (i, &thickness) in thicknesses.iter().enumerate() {
    reflectivity_spectra.row_mut(i).assign(&transfer_matrix_reflectivity(n_air, &n_thin_film, n_si, thickness, &wavelengths));
}

// Plot reflectivity spectra for selected thicknesses
let root = BitMapBackend::new("reflectivity_spectra.png", (1200, 800)).into_drawing_area();
root.fill(&WHITE).unwrap();
let mut chart = ChartBuilder::on(&root)
    .caption("Reflectivity Spectra of stack", ("sans-serif", 50).into_font())
    .margin(10)
    .x_label_area_size(30)
    .y_label_area_size(30)
    .build_cartesian_2d(200.0..800.0, 0.0..1.0)
    .unwrap();

chart.configure_mesh().draw().unwrap();

for (i, &thickness) in thicknesses.iter().step_by(1000).enumerate() {
    chart
        .draw_series(LineSeries::new(
            wavelengths.iter().zip(reflectivity_spectra.row(i).iter()).map(|(&x, &y)| (x, y)),
            &Palette99::pick(i),
        ))
        .unwrap()
        .label(format!("Thickness = {} Å", thickness))
        .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &Palette99::pick(i)));
}

chart.configure_series_labels().background_style(&WHITE.mix(0.8)).border_style(&BLACK).draw().unwrap();