'use strict';
import * as arrayStat from './array';
/**
 *
 * @param {number} a
 * @param {number} b
 * @returns
 */
function compareNumbers(a, b) {
    return a - b;
}
/**
 *
 * @param {number[][]} matrix
 * @returns
 */
export function max(matrix) {
    let max = -Infinity;
    for (let i = 0; i < matrix.length; i++) {
        for (let j = 0; j < matrix[i].length; j++) {
            if (matrix[i][j] > max)
                max = matrix[i][j];
        }
    }
    return max;
}
/**
 *
 * @param {number[][]} matrix
 * @returns
 */
export function min(matrix) {
    let min = Infinity;
    for (let i = 0; i < matrix.length; i++) {
        for (let j = 0; j < matrix[i].length; j++) {
            if (matrix[i][j] < min)
                min = matrix[i][j];
        }
    }
    return min;
}
/**
 *
 * @param {number[][]} matrix
 */
export function minMax(matrix) {
    let min = Infinity;
    let max = -Infinity;
    for (let i = 0; i < matrix.length; i++) {
        for (let j = 0; j < matrix[i].length; j++) {
            if (matrix[i][j] < min)
                min = matrix[i][j];
            if (matrix[i][j] > max)
                max = matrix[i][j];
        }
    }
    return {
        min: min,
        max: max
    };
}
/**
 *
 * @param {number[][]} matrix
 * @param eps = 0
 * @returns
 */
export function entropy(matrix, eps = 0) {
    let sum = 0, l1 = matrix.length, l2 = matrix[0].length;
    for (let i = 0; i < l1; i++) {
        for (let j = 0; j < l2; j++) {
            sum += matrix[i][j] * Math.log(matrix[i][j] + eps);
        }
    }
    return -sum;
}
/**
 *
 * @param {number[][]} matrix
 * @param dimension = 0
 * @returns
 */
export function mean(matrix, dimension = 0) {
    let rows = matrix.length, cols = matrix[0].length, theMean, N, i, j;
    if (dimension === -1) {
        theMean = [0];
        N = rows * cols;
        for (i = 0; i < rows; i++) {
            for (j = 0; j < cols; j++) {
                theMean[0] += matrix[i][j];
            }
        }
        theMean[0] /= N;
    }
    else if (dimension === 0) {
        theMean = new Array(cols);
        N = rows;
        for (j = 0; j < cols; j++) {
            theMean[j] = 0;
            for (i = 0; i < rows; i++) {
                theMean[j] += matrix[i][j];
            }
            theMean[j] /= N;
        }
    }
    else if (dimension === 1) {
        theMean = new Array(rows);
        N = cols;
        for (j = 0; j < rows; j++) {
            theMean[j] = 0;
            for (i = 0; i < cols; i++) {
                theMean[j] += matrix[j][i];
            }
            theMean[j] /= N;
        }
    }
    else {
        throw new Error('Invalid dimension');
    }
    return theMean;
}
/**
 *
 * @param {number[][]} matrix
 * @param dimension = 0
 * @returns
 */
export function sum(matrix, dimension = 0) {
    let rows = matrix.length, cols = matrix[0].length, theSum;
    if (dimension === -1) {
        theSum = [0];
        for (let i = 0; i < rows; i++) {
            for (let j = 0; j < cols; j++) {
                theSum[0] += matrix[i][j];
            }
        }
    }
    else if (dimension === 0) {
        theSum = new Array(cols);
        for (let j = 0; j < cols; j++) {
            theSum[j] = 0;
            for (let i = 0; i < rows; i++) {
                theSum[j] += matrix[i][j];
            }
        }
    }
    else if (dimension === 1) {
        theSum = new Array(rows);
        for (let j = 0; j < rows; j++) {
            theSum[j] = 0;
            for (let i = 0; i < cols; i++) {
                theSum[j] += matrix[j][i];
            }
        }
    }
    else {
        throw new Error('Invalid dimension');
    }
    return theSum;
}
/**
 *
 * @param {number[][]} matrix
 * @param dimension = 0
 * @returns
 */
export function product(matrix, dimension = 0) {
    let rows = matrix.length, cols = matrix[0].length, theProduct;
    if (dimension === -1) {
        theProduct = [1];
        for (let i = 0; i < rows; i++) {
            for (let j = 0; j < cols; j++) {
                theProduct[0] *= matrix[i][j];
            }
        }
    }
    else if (dimension === 0) {
        theProduct = new Array(cols);
        for (let j = 0; j < cols; j++) {
            theProduct[j] = 1;
            for (let i = 0; i < rows; i++) {
                theProduct[j] *= matrix[i][j];
            }
        }
    }
    else if (dimension === 1) {
        theProduct = new Array(rows);
        for (let j = 0; j < rows; j++) {
            theProduct[j] = 1;
            for (let i = 0; i < cols; i++) {
                theProduct[j] *= matrix[j][i];
            }
        }
    }
    else {
        throw new Error('Invalid dimension');
    }
    return theProduct;
}
/**
 *
 * @param {number[][]} matrix
 * @param {number[]} means?
 * @param {boolean} unbiased?
 * @returns
 */
export function standardDeviation(matrix, unbiased, means) {
    let vari = variance(matrix, means, unbiased), l = vari.length;
    for (let i = 0; i < l; i++) {
        vari[i] = Math.sqrt(vari[i]);
    }
    return vari;
}
/**
 *
 * @param {number[][]} matrix
 * @param {number[]} means
 * @param unbiased = true
 * @returns
 */
export function variance(matrix, means = mean(matrix), unbiased = true) {
    let rows = matrix.length;
    if (rows === 0)
        return [];
    let cols = matrix[0].length;
    let vari = new Array(cols);
    for (let j = 0; j < cols; j++) {
        let sum1 = 0, sum2 = 0;
        //x = 0;
        for (let i = 0; i < rows; i++) {
            let x = matrix[i][j] - means[j];
            sum1 += x;
            sum2 += x * x;
        }
        if (unbiased) {
            vari[j] = (sum2 - ((sum1 * sum1) / rows)) / (rows - 1);
        }
        else {
            vari[j] = (sum2 - ((sum1 * sum1) / rows)) / rows;
        }
    }
    return vari;
}
/**
 *
 * @param {number[][]} matrix
 * @returns
 */
export function median(matrix) {
    let rows = matrix.length, cols = matrix[0].length;
    let medians = new Array(cols);
    for (let i = 0; i < cols; i++) {
        let data = new Array(rows);
        for (let j = 0; j < rows; j++) {
            data[j] = matrix[j][i];
        }
        data.sort(compareNumbers);
        let N = data.length;
        if (N % 2 === 0) {
            medians[i] = (data[N / 2] + data[(N / 2) - 1]) * 0.5;
        }
        else {
            medians[i] = data[Math.floor(N / 2)];
        }
    }
    return medians;
}
/**
 *
 * @param {number[][]} matrix
 * @returns
 */
export function mode(matrix) {
    let rows = matrix.length, cols = matrix[0].length, modes = new Array(cols);
    for (let i = 0; i < cols; i++) {
        let itemCount = new Array(rows);
        for (let k = 0; k < rows; k++) {
            itemCount[k] = 0;
        }
        let itemArray = new Array(rows);
        let count = 0;
        for (let j = 0; j < rows; j++) {
            let index = itemArray.indexOf(matrix[j][i]);
            if (index >= 0) {
                itemCount[index]++;
            }
            else {
                itemArray[count] = matrix[j][i];
                itemCount[count] = 1;
                count++;
            }
        }
        let maxValue = 0, maxIndex = 0;
        for (let j = 0; j < count; j++) {
            if (itemCount[j] > maxValue) {
                maxValue = itemCount[j];
                maxIndex = j;
            }
        }
        modes[i] = itemArray[maxIndex];
    }
    return modes;
}
/**
 *
 * @param matrix
 * @param unbiased
 */
export function skewness(matrix, unbiased = true) {
    let means = mean(matrix);
    let n = matrix.length, l = means.length;
    let skew = new Array(l);
    for (let j = 0; j < l; j++) {
        let s2 = 0, s3 = 0;
        for (let i = 0; i < n; i++) {
            let dev = matrix[i][j] - means[j];
            s2 += dev * dev;
            s3 += dev * dev * dev;
        }
        let m2 = s2 / n;
        let m3 = s3 / n;
        let g = m3 / Math.pow(m2, 3 / 2);
        if (unbiased) {
            let a = Math.sqrt(n * (n - 1));
            let b = n - 2;
            skew[j] = (a / b) * g;
        }
        else {
            skew[j] = g;
        }
    }
    return skew;
}
/**
 *
 * @param matrix
 * @param unbiased
 */
export function kurtosis(matrix, unbiased = true) {
    let means = mean(matrix);
    let n = matrix.length, m = matrix[0].length;
    let kurt = new Array(m);
    for (let j = 0; j < m; j++) {
        let s2 = 0, s4 = 0;
        for (let i = 0; i < n; i++) {
            let dev = matrix[i][j] - means[j];
            s2 += dev * dev;
            s4 += dev * dev * dev * dev;
        }
        let m2 = s2 / n;
        let m4 = s4 / n;
        if (unbiased) {
            let v = s2 / (n - 1);
            let a = (n * (n + 1)) / ((n - 1) * (n - 2) * (n - 3));
            let b = s4 / (v * v);
            let c = ((n - 1) * (n - 1)) / ((n - 2) * (n - 3));
            kurt[j] = a * b - 3 * c;
        }
        else {
            kurt[j] = m4 / (m2 * m2) - 3;
        }
    }
    return kurt;
}
/**
 *
 * @param matrix
 */
export function standardError(matrix) {
    let samples = matrix.length;
    let standardDeviations = standardDeviation(matrix);
    let l = standardDeviations.length;
    let standardErrors = new Array(l);
    let sqrtN = Math.sqrt(samples);
    for (let i = 0; i < l; i++) {
        standardErrors[i] = standardDeviations[i] / sqrtN;
    }
    return standardErrors;
}
/**
 *
 * @param {number[][]} matrix
 * @param {number} dimension?
 * @returns
 */
export function covariance(matrix, dimension) {
    return scatter(matrix, undefined, dimension);
}
/**
 *
 * @param {number[][]} matrix
 * @param {number} divisor?
 * @param dimension = 0
 * @returns
 */
export function scatter(matrix, divisor, dimension = 0) {
    if (typeof (divisor) === 'undefined') {
        if (dimension === 0) {
            divisor = matrix.length - 1;
        }
        else if (dimension === 1) {
            divisor = matrix[0].length - 1;
        }
    }
    let means = mean(matrix, dimension);
    let rows = matrix.length;
    if (rows === 0) {
        return [[]];
    }
    let cols = matrix[0].length, cov, j, s;
    if (dimension === 0) {
        cov = new Array(cols);
        for (let i = 0; i < cols; i++) {
            cov[i] = new Array(cols);
        }
        for (let i = 0; i < cols; i++) {
            for (j = i; j < cols; j++) {
                s = 0;
                for (let k = 0; k < rows; k++) {
                    s += (matrix[k][j] - means[j]) * (matrix[k][i] - means[i]);
                }
                s /= divisor;
                cov[i][j] = s;
                cov[j][i] = s;
            }
        }
    }
    else if (dimension === 1) {
        cov = new Array(rows);
        for (let i = 0; i < rows; i++) {
            cov[i] = new Array(rows);
        }
        for (let i = 0; i < rows; i++) {
            for (j = i; j < rows; j++) {
                s = 0;
                for (let k = 0; k < cols; k++) {
                    s += (matrix[j][k] - means[j]) * (matrix[i][k] - means[i]);
                }
                s /= divisor;
                cov[i][j] = s;
                cov[j][i] = s;
            }
        }
    }
    else {
        throw new Error('Invalid dimension');
    }
    return cov;
}
/**
 *
 * @param {number[][]} matrix
 * @returns
 */
export function correlation(matrix) {
    let means = mean(matrix), standardDeviations = standardDeviation(matrix, true, means), scores = zScores(matrix, means, standardDeviations), rows = matrix.length, cols = matrix[0].length;
    let cor = new Array(cols);
    for (let i = 0; i < cols; i++) {
        cor[i] = new Array(cols);
    }
    for (let i = 0; i < cols; i++) {
        for (let j = i; j < cols; j++) {
            let c = 0;
            for (let k = 0, l = scores.length; k < l; k++) {
                c += scores[k][j] * scores[k][i];
            }
            c /= rows - 1;
            cor[i][j] = c;
            cor[j][i] = c;
        }
    }
    return cor;
}
/**
 *
 * @param matrix
 * @param means
 * @param standardDeviations
 */
export function zScores(matrix, means = mean(matrix), standardDeviations = standardDeviation(matrix, true, means)) {
    return standardize(center(matrix, means, false), standardDeviations, true);
}
/**
 *
 * @param {number[][]} matrix
 * @param {number[]} means
 * @param {boolean} inPlace
 * @returns
 */
export function center(matrix, means = mean(matrix), inPlace) {
    let result = matrix, l = matrix.length;
    if (!inPlace) {
        result = new Array(l);
        for (let i = 0; i < l; i++) {
            result[i] = new Array(matrix[i].length);
        }
    }
    for (let i = 0; i < l; i++) {
        let row = result[i];
        for (let j = 0, jj = row.length; j < jj; j++) {
            row[j] = matrix[i][j] - means[j];
        }
    }
    return result;
}
/**
 *
 * @param {number[][]} matrix
 * @param standardDeviations = standardDeviation(matrix)
 * @param {boolean} inPlace
 * @returns
 */
export function standardize(matrix, standardDeviations = standardDeviation(matrix), inPlace) {
    let result = matrix, l = matrix.length;
    if (!inPlace) {
        result = new Array(l);
        for (let i = 0; i < l; i++) {
            result[i] = new Array(matrix[i].length);
        }
    }
    for (let i = 0; i < l; i++) {
        let resultRow = result[i];
        let sourceRow = matrix[i];
        for (let j = 0, jj = resultRow.length; j < jj; j++) {
            if (standardDeviations[j] !== 0 && !isNaN(standardDeviations[j])) {
                resultRow[j] = sourceRow[j] / standardDeviations[j];
            }
        }
    }
    return result;
}
/**
 *
 * @param {number[][]} matrix
 * @param {number[]} weights
 * @returns
 */
export function weightedVariance(matrix, weights) {
    let means = mean(matrix);
    let rows = matrix.length;
    if (rows === 0)
        return [];
    let cols = matrix[0].length;
    let vari = new Array(cols);
    for (let j = 0; j < cols; j++) {
        let sum = 0;
        let a = 0, b = 0;
        for (let i = 0; i < rows; i++) {
            let z = matrix[i][j] - means[j];
            let w = weights[i];
            sum += w * (z * z);
            b += w;
            a += w * w;
        }
        vari[j] = sum * (b / (b * b - a));
    }
    return vari;
}
/**
 *
 * @param {number[][]} matrix
 * @param {number[]} weights
 * @param dimension = 0
 * @returns
 */
export function weightedMean(matrix, weights, dimension = 0) {
    let rows = matrix.length;
    if (rows === 0) {
        return [];
    }
    let cols = matrix[0].length, means, w, row;
    if (dimension === 0) {
        means = new Array(cols);
        for (let i = 0; i < cols; i++) {
            means[i] = 0;
        }
        for (let i = 0; i < rows; i++) {
            row = matrix[i];
            w = weights[i];
            for (let j = 0; j < cols; j++) {
                means[j] += row[j] * w;
            }
        }
    }
    else if (dimension === 1) {
        means = new Array(rows);
        for (let i = 0; i < rows; i++) {
            means[i] = 0;
        }
        for (let j = 0; j < rows; j++) {
            row = matrix[j];
            w = weights[j];
            for (let i = 0; i < cols; i++) {
                means[j] += row[i] * w;
            }
        }
    }
    else {
        throw new Error('Invalid dimension');
    }
    let weightSum = arrayStat.sum(weights);
    if (weightSum !== 0) {
        for (let i = 0, ii = means.length; i < ii; i++) {
            means[i] /= weightSum;
        }
    }
    return means;
}
/**
 *
 * @param {number[][]} matrix
 * @param {number[]} weights
 * @param dimension = 0
 * @param means = weightedMean(matrix
 * @param weights
 * @param dimension)
 * @returns
 */
export function weightedCovariance(matrix, weights, dimension = 0, means = weightedMean(matrix, weights, dimension)) {
    let s1 = 0, s2 = 0;
    for (let i = 0, ii = weights.length; i < ii; i++) {
        s1 += weights[i];
        s2 += weights[i] * weights[i];
    }
    let factor = s1 / (s1 * s1 - s2);
    return weightedScatter(matrix, weights, means, factor, dimension);
}
/**
 *
 * @param {number[][]} matrix
 * @param {number[]} weights
 * @param {number[]} means
 * @param factor = 1
 * @param {number} dimension?
 * @returns
 */
export function weightedScatter(matrix, weights, means, factor = 1, dimension = 0) {
    means = means || weightedMean(matrix, weights, dimension);
    let rows = matrix.length;
    if (rows === 0) {
        return [[]];
    }
    let cols = matrix[0].length, cov;
    if (dimension === 0) {
        cov = new Array(cols);
        for (let i = 0; i < cols; i++) {
            cov[i] = new Array(cols);
        }
        for (let i = 0; i < cols; i++) {
            for (let j = i; j < cols; j++) {
                let s = 0;
                for (let k = 0; k < rows; k++) {
                    s += weights[k] * (matrix[k][j] - means[j]) * (matrix[k][i] - means[i]);
                }
                cov[i][j] = s * factor;
                cov[j][i] = s * factor;
            }
        }
    }
    else if (dimension === 1) {
        cov = new Array(rows);
        for (let i = 0; i < rows; i++) {
            cov[i] = new Array(rows);
        }
        for (let i = 0; i < rows; i++) {
            for (let j = i; j < rows; j++) {
                let s = 0;
                for (let k = 0; k < cols; k++) {
                    s += weights[k] * (matrix[j][k] - means[j]) * (matrix[i][k] - means[i]);
                }
                cov[i][j] = s * factor;
                cov[j][i] = s * factor;
            }
        }
    }
    else {
        throw new Error('Invalid dimension');
    }
    return cov;
}
//# sourceMappingURL=matrix.js.map