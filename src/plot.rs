use crate::utils;
use plotters::prelude::*;

/// Boundary box delineating extent of plotting data.
pub struct Bbox {
    xmin: f64,
    xmax: f64,
    ymin: f64,
    ymax: f64,
}

impl Bbox {
    fn new(x: &[f64], y: &[f64]) -> Self {
        let xmin: f64 = x.iter().cloned().fold(0.0, f64::min);
        let xmax: f64 = x.iter().cloned().fold(f64::NAN, f64::max);
        let ymin: f64 = y.iter().cloned().fold(0.0, f64::min);
        let ymax: f64 = y.iter().cloned().fold(f64::NAN, f64::max);
        Bbox {
            xmin,
            xmax,
            ymin,
            ymax
        }
    }
}

/// Generic function for plotting results on the fly.
pub fn xy(x: &[f64], y: &[f64], path: &str) -> Result<(), Box<dyn std::error::Error>> {
    let bbox = Bbox::new(&x, &y);
    let xy: Vec<(f64, f64)> = x.iter().cloned().zip(y.iter().cloned()).collect();
    let root = BitMapBackend::new(path, (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    root.margin(10, 10, 10, 10);
    let mut chart = ChartBuilder::on(&root)
               // .caption("Title", ("sans-serif", 16).into_font())
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(bbox.xmin..bbox.xmax, bbox.ymin..bbox.ymax)?;
    chart
        .configure_mesh()
        // .x_labels(5)
        // .y_labels(5)
        .y_label_formatter(&|x| format!("{:.2}", x))
        .x_label_formatter(&|x| format!("{:.2}", x))
        // .x_desc("Value")
        // .y_desc("CDF")
        .draw()?;
    chart.draw_series(xy.iter().map(|x| Circle::new((x.0, x.1), 2, BLUE.filled())))?;

    Ok(())


}

/// Compare the CDF of two accumulation records.
pub fn comp_cdf(a: &[f64], b: &[f64], title: &str) -> Result<(), Box<dyn std::error::Error>> {
    let a = utils::cdf(a);
    let b = utils::cdf(b);
    let mut ab = a.clone();
    let mut c = b.clone();
    ab.append(&mut c);

    let ymin = ab.iter().map(|xi| xi.1).fold(0.0, f64::min);
    let ymax = ab.iter().map(|xi| xi.1).fold(0.0, f64::max);
    let xmin = ab.iter().map(|xi| xi.0).fold(0.0, f64::min);
    let xmax = ab.iter().map(|xi| xi.0).fold(0.0, f64::max);
    let root = BitMapBackend::new(title, (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    root.margin(10, 10, 10, 10);
    // construct a chart context
    let mut chart = ChartBuilder::on(&root)
        // Set the caption of the chart
        //        .caption("Title", ("sans-serif", 40).into_font())
        // Set the size of the label region
        .x_label_area_size(40)
        .y_label_area_size(60)
        // Finally attach a coordinate on the drawing area and make a chart context
        .build_cartesian_2d(xmin..xmax, ymin..ymax)?;

    // Then we can draw a mesh
    chart
        .configure_mesh()
        // We can customize the maximum number of labels allowed for each axis
        .x_labels(5)
        .y_labels(5)
        // We can also change the format of the label text
        .y_label_formatter(&|x| format!("{:.2}", x))
        .x_label_formatter(&|x| format!("{:.2}", x))
        .x_desc("Value")
        .y_desc("CDF")
        .draw()?;

    // And we can draw something in the drawing area
    chart
        .draw_series(LineSeries::new(a.clone(), &BLACK))?
        .label("synthetic")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));
    // Similarly, we can draw point series
    chart.draw_series(PointSeries::of_element(a, 2, &BLUE, &|c, s, st| {
        return EmptyElement::at(c)    // We want to construct a composed element on-the-fly
                + Circle::new((0, 0), s, st.filled()); // At this point, the new pixel coordinate is established
                                                       //                + Text::new(format!("{:?}", c), (10, 0), ("sans-serif", 10).into_font());
    }))?;

    chart
        .draw_series(LineSeries::new(b.clone(), &BLACK))?
        .label("observed")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));
    chart.draw_series(PointSeries::of_element(b, 2, &GREEN, &|c, s, st| {
        return EmptyElement::at(c)    // We want to construct a composed element on-the-fly
                + Circle::new((0, 0), s, st.filled()); // At this point, the new pixel coordinate is established
                                                       //                + Text::new(format!("{:?}", c), (10, 0), ("sans-serif", 10).into_font());
    }))?;
    chart
        .configure_series_labels()
        .background_style(WHITE.filled())
        .draw()?;
    Ok(())
}

/// Scatter plot for comparing goodness-of-fit statistics from model runs.
pub fn scatter(x: &[f64], y: &[f64], title: &str) -> Result<(), Box<dyn std::error::Error>> {
    let xy: Vec<(f64, f64)> = x.iter().cloned().zip(y.iter().cloned()).collect();

    let ymin = xy.iter().map(|xi| xi.1).fold(1.0, f64::min);
    let ymax = xy.iter().map(|xi| xi.1).fold(0.0, f64::max);
    let xmin = xy.iter().map(|xi| xi.0).fold(1.0, f64::min);
    let xmax = xy.iter().map(|xi| xi.0).fold(0.0, f64::max);
    let root = BitMapBackend::new(title, (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    root.margin(10, 10, 10, 10);
    let mut chart = ChartBuilder::on(&root)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(xmin..xmax, ymin..ymax)?;

    chart
        .configure_mesh()
        .x_labels(5)
        .y_labels(5)
        .y_label_formatter(&|x| format!("{:.2}", x))
        .x_label_formatter(&|x| format!("{:.2}", x))
        .x_desc("Rate")
        .y_desc("Fit")
        .draw()?;
    chart.draw_series(xy.iter().map(|x| Circle::new((x.0, x.1), 2, BLUE.filled())))?;

    Ok(())
}

/// Experimental box-and-whisker plot of transit times.
pub fn whisker_for_facies(
    df: &plotters::data::Quartiles,
    ff: &plotters::data::Quartiles,
    fg: &plotters::data::Quartiles,
    title: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let values_range = plotters::data::fitting_range(
        df.values()
            .iter()
            .chain(ff.values().iter().chain(fg.values().iter())),
    );

    let facies_axis = ["debris flows", "fines", "gravels"];
    let root = BitMapBackend::new(title, (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    root.margin(10, 10, 10, 10);
    // construct a chart context
    let mut chart = ChartBuilder::on(&root)
        // Set the caption of the chart
        //        .caption("Title", ("sans-serif", 40).into_font())
        // Set the size of the label region
        .x_label_area_size(40)
        .y_label_area_size(60)
        // Finally attach a coordinate on the drawing area and make a chart context
        .build_cartesian_2d(facies_axis[..].into_segmented(), values_range)?;

    // Then we can draw a mesh
    chart
        .configure_mesh()
        // We can customize the maximum number of labels allowed for each axis
        .y_labels(5)
        // We can also change the format of the label text
        .y_label_formatter(&|x| format!("{:.2}", x))
        .x_desc("Facies")
        .y_desc("Transit Time")
        .draw()?;

    // And we can draw something in the drawing area

    chart.draw_series(vec![
        Boxplot::new_vertical(SegmentValue::CenterOf(&"debris flows"), &df),
        Boxplot::new_vertical(SegmentValue::CenterOf(&"fines"), &ff),
        Boxplot::new_vertical(SegmentValue::CenterOf(&"gravels"), &fg),
    ])?;

    chart
        .configure_series_labels()
        .background_style(WHITE.filled())
        .draw()?;
    Ok(())
}
