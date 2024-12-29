pub fn eval(message: &str) -> Vec<u8> {
    message.as_bytes().to_vec()
}
pub fn erasure_coding(message: &str) -> Vec<(f64, f64)> {
    let message_bytes = eval(message);
    let polynomial_coeffs = message_bytes
        .clone()
        .into_iter()
        .map(|f| f as f64)
        .collect::<Vec<f64>>();
    let mut code_word: Vec<f64> = vec![];
    for i in 0..polynomial_coeffs.len() {
        code_word.push(poly_eval(&polynomial_coeffs, (i + 1) as f64));
    }
    let parity_1 = poly_eval(&polynomial_coeffs, (message_bytes.clone().len() + 1) as f64);
    let parity_2 = poly_eval(&polynomial_coeffs, (message_bytes.clone().len() + 2) as f64);
    let parity_3 = poly_eval(&polynomial_coeffs, (message_bytes.len() + 3) as f64);

    code_word.push(parity_1);
    code_word.push(parity_2);
    code_word.push(parity_3);

    let mut result: Vec<(f64, f64)> = vec![];
    for i in 0..code_word.len() {
        result.push(((i + 1) as f64, code_word[i]));
    }
    result
}

pub fn lagrange_interpolation(points: Vec<(f64,f64)>) -> Vec<f64> {
    // let mut points = vec![];
    // for i in 0..message.len() {
    //     points.push(((i + 1) as f64, message[i]));
    // }
    println!("{:?}", points);
    let n = points.len();
    let mut coeffs = vec![0.0; n]; // Resultant polynomial coefficients

    for i in 0..n {
        let (x_i, y_i) = points[i];
        let mut basis_coeffs = vec![1.0]; // Coefficients for L_i(x)

        for j in 0..n {
            if i != j {
                let (x_j, _) = points[j];

                // Update basis polynomial for L_i(x)
                let mut new_basis = vec![0.0; basis_coeffs.len() + 1];
                for k in 0..basis_coeffs.len() {
                    new_basis[k] -= basis_coeffs[k] * x_j; // -x_j
                    new_basis[k + 1] += basis_coeffs[k]; // +1 * x
                }
                basis_coeffs = new_basis;

                // Scale by the denominator (x_i - x_j)
                let scale = 1.0 / (x_i - x_j);
                for k in 0..basis_coeffs.len() {
                    basis_coeffs[k] *= scale;
                }
            }
        }

        // Add y_i * L_i(x) to the result
        for k in 0..basis_coeffs.len() {
            coeffs[k] += y_i * basis_coeffs[k];
        }
    }

    coeffs
}

pub fn poly_eval(coeffs: &Vec<f64>, point: f64) -> f64 {
    let mut result = 0.0;

    for i in 0..coeffs.len() {
        result += coeffs[i] * point.powi(i.try_into().unwrap());
    }
    result
}

pub fn recover(coeffs: &Vec<(f64,f64)>,point: f64)->f64{
    let coeff_y: Vec<_> = coeffs.into_iter().map(|v| v.1).collect();

    let mut result = 0.0;

    for i in 0..coeff_y.len() {
        result += coeff_y[i] * point.powi(i.try_into().unwrap());
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let result = erasure_coding("Hello");

        println!("{:?}", result);
    }

    // #[test]
    // fn lagrange_works() {
    //     let result = lagrange_interpolation(vec![72.0, 101.0, 108.0, 108.0, 111.0]);
    //     let d = result.len();

    //     // println!("{:?}", result);
    //     assert_eq!(result[0], 1.0);
    //     assert_eq!(result[1], 99.91666666666654);
    //     assert_eq!(result[2], -33.291666666666515);
    //     assert_eq!(result[3], 4.5833333333333215);
    //     assert_eq!(result[4], -0.20833333333333215);

    //     let eval_1 = poly_eval(&result.clone(), (d + 1) as f64);
    //     let eval_2 = poly_eval(&result.clone(), (d + 2) as f64);
    //     let eval_3 = poly_eval(&result, (d + 3) as f64);

    //     assert_eq!(eval_1, 122.00000000000381);
    //     assert_eq!(eval_2, 141.0000000000053);
    //     assert_eq!(eval_3, 163.0000000000075);
    // }
    
    
    #[test]
    fn erasure_coding_works() {
        let result = erasure_coding("Hello");
        let malformed_result = vec![(1.0, 500.0), (3.0, 13254.0), (5.0, 86152.0), (6.0, 171750.0), (7.0, 309626.0)];
        let recovered_polynomial=lagrange_interpolation(malformed_result);

        println!("Malformed:{:?}", recovered_polynomial);

        println!("At 2:{}",poly_eval(&recovered_polynomial, 2.0));
    }
}
