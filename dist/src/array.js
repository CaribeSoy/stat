'use strict';
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
 * Computes the sum of the given values
 * @param {Array} values
 * @returns {number}
 */
export function sum(values) {
    let total = 0;
    for (let i = 0; i < values.length; i++) {
        total += values[i];
    }
    return total;
}
/**
 * Computes the maximum of the given values
 * @param {Array} values
 * @returns {number}
 */
export function max(values) {
    let max = values[0];
    const l = values.length;
    for (let i = 1; i < l; i++) {
        if (values[i] > max) {
            max = values[i];
        }
    }
    return max;
}
/**
 * Computes the minimum of the given values
 * @param {Array} values
 * @returns {number}
 */
export function min(values) {
    let min = values[0];
    let l = values.length;
    for (let i = 1; i < l; i++) {
        if (values[i] < min) {
            min = values[i];
        }
    }
    return min;
}
/**
 * Computes the min and max of the given values
 * @param {Array} values
 * @returns {{min: number, max: number}}
 */
export function minMax(values) {
    let min = values[0];
    let max = values[0];
    const l = values.length;
    for (let i = 1; i < l; i++) {
        if (values[i] < min)
            min = values[i];
        if (values[i] > max)
            max = values[i];
    }
    return {
        min: min,
        max: max
    };
}
/**
 * Computes the arithmetic mean of the given values
 * @param {Array} values
 * @returns {number}
 */
export function arithmeticMean(values) {
    let sum = 0;
    const l = values.length;
    for (let i = 0; i < l; i++) {
        sum += values[i];
    }
    return sum / l;
}
/**
 * {@link arithmeticMean}
 */
export var mean = arithmeticMean;
/**
 * Computes the geometric mean of the given values
 * @param {Array} values
 * @returns {number}
 */
export function geometricMean(values) {
    let mul = 1;
    const l = values.length;
    for (let i = 0; i < l; i++) {
        mul *= values[i];
    }
    return Math.pow(mul, 1 / l);
}
/**
 * Computes the mean of the log of the given values
 * If the return value is exponentiated, it gives the same result as the
 * geometric mean.
 * @param {Array} values
 * @returns {number}
 */
export function logMean(values) {
    let lnsum = 0;
    const l = values.length;
    for (let i = 0; i < l; i++) {
        lnsum += Math.log(values[i]);
    }
    return lnsum / l;
}
/**
 * Computes the weighted grand mean for a list of means and sample sizes
 * @param {Array} means - Mean values for each set of samples
 * @param {Array} samples - Number of original values for each set of samples
 * @returns {number}
 */
export function grandMean(means, samples) {
    let sum = 0;
    let n = 0;
    const l = means.length;
    for (let i = 0; i < l; i++) {
        sum += samples[i] * means[i];
        n += samples[i];
    }
    return sum / n;
}
/**
 * Computes the truncated mean of the given values using a given percentage
 * @param {Array} values
 * @param {number} percent - The percentage of values to keep (range: [0,1])
 * @param {boolean} [alreadySorted=false]
 * @returns {number}
 */
export function truncatedMean(values, percent, alreadySorted = false) {
    if (!alreadySorted) {
        values = [].concat(values).sort(compareNumbers);
    }
    const l = values.length;
    const k = Math.floor(l * percent);
    let sum = 0;
    for (let i = k; i < (l - k); i++) {
        sum += values[i];
    }
    return sum / (l - 2 * k);
}
/**
 * Computes the harmonic mean of the given values
 * @param {Array} values
 * @returns {number}
 */
export function harmonicMean(values) {
    let sum = 0;
    const l = values.length;
    for (let i = 0; i < l; i++) {
        if (values[i] === 0) {
            throw new RangeError('value at index ' + i + 'is zero');
        }
        sum += 1 / values[i];
    }
    return l / sum;
}
/**
 * Computes the contra-harmonic mean of the given values
 * @param {Array} values
 * @returns {number}
 */
export function contraHarmonicMean(values) {
    let r1 = 0;
    let r2 = 0;
    const l = values.length;
    for (let i = 0; i < l; i++) {
        r1 += values[i] * values[i];
        r2 += values[i];
    }
    if (r2 < 0) {
        throw new RangeError('sum of values is negative');
    }
    return r1 / r2;
}
/**
 * Computes the median of the given values
 * @param {Array} values
 * @param {boolean} [alreadySorted=false]
 * @returns {number}
 */
export function median(values, alreadySorted = false) {
    if (!alreadySorted) {
        values = [].concat(values).sort(compareNumbers);
    }
    const l = values.length;
    const half = Math.floor(l / 2);
    if (l % 2 === 0) {
        return (values[half - 1] + values[half]) * 0.5;
    }
    else {
        return values[half];
    }
}
/**
 * Computes the variance of the given values
 * @param {Array} values
 * @param {boolean} [unbiased=true] - if true, divide by (n-1); if false, divide by n.
 * @returns {number}
 */
export function variance(values, unbiased = true) {
    const theMean = mean(values);
    let theVariance = 0;
    const l = values.length;
    for (let i = 0; i < l; i++) {
        var x = values[i] - theMean;
        theVariance += x * x;
    }
    if (unbiased) {
        return theVariance / (l - 1);
    }
    else {
        return theVariance / l;
    }
}
/**
 * Computes the standard deviation of the given values
 * @param {Array} values
 * @param {boolean} [unbiased=true] - if true, divide by (n-1); if false, divide by n.
 * @returns {number}
 */
export function standardDeviation(values, unbiased = true) {
    return Math.sqrt(variance(values, unbiased));
}
export function standardError(values) {
    return standardDeviation(values) / Math.sqrt(values.length);
}
/**
 * IEEE Transactions on biomedical engineering, vol. 52, no. 1, January 2005, p. 76-
 * Calculate the standard deviation via the Median of the absolute deviation
 *  The formula for the standard deviation only holds for Gaussian random variables.
 * @returns {{mean: number, stdev: number}}
 */
export function robustMeanAndStdev(y) {
    let mean = 0, stdev;
    const length = y.length;
    for (let i = 0; i < length; i++) {
        mean += y[i];
    }
    mean /= length;
    let averageDeviations = new Array(length);
    for (let i = 0; i < length; i++) {
        averageDeviations[i] = Math.abs(y[i] - mean);
    }
    averageDeviations.sort(compareNumbers);
    if (length % 2 === 1) {
        stdev = averageDeviations[(length - 1) / 2] / 0.6745;
    }
    else {
        stdev = 0.5 * (averageDeviations[length / 2] + averageDeviations[length / 2 - 1]) / 0.6745;
    }
    return {
        mean: mean,
        stdev: stdev
    };
}
export function quartiles(values, alreadySorted = false) {
    if (!alreadySorted) {
        values = [].concat(values).sort(compareNumbers);
    }
    const quart = values.length / 4;
    const q1 = values[Math.ceil(quart) - 1];
    const q2 = median(values, true);
    const q3 = values[Math.ceil(quart * 3) - 1];
    return { q1: q1, q2: q2, q3: q3 };
}
export function pooledStandardDeviation(samples, unbiased = true) {
    return Math.sqrt(pooledVariance(samples, unbiased));
}
export function pooledVariance(samples, unbiased = true) {
    let sum = 0;
    let length = 0;
    const l = samples.length;
    for (let i = 0; i < l; i++) {
        let values = samples[i];
        let vari = variance(values);
        sum += (values.length - 1) * vari;
        if (unbiased) {
            length += values.length - 1;
        }
        else {
            length += values.length;
        }
    }
    return sum / length;
}
/**
 *
 * @param {number[]} values
 * @returns
 */
export function mode(values) {
    const l = values.length;
    let itemCount = new Array(l);
    for (let i = 0; i < l; i++) {
        itemCount[i] = 0;
    }
    let itemArray = new Array(l);
    let count = 0;
    for (let i = 0; i < l; i++) {
        let index = itemArray.indexOf(values[i]);
        if (index >= 0)
            itemCount[index]++;
        else {
            itemArray[count] = values[i];
            itemCount[count] = 1;
            count++;
        }
    }
    let maxValue = 0, maxIndex = 0;
    for (let i = 0; i < count; i++) {
        if (itemCount[i] > maxValue) {
            maxValue = itemCount[i];
            maxIndex = i;
        }
    }
    return itemArray[maxIndex];
}
/**
 *
 * @param {number[]} vector1
 * @param {number[]} vector2
 * @param unbiased = true
 * @returns
 */
export function covariance(vector1, vector2, unbiased = true) {
    let mean1 = mean(vector1);
    let mean2 = mean(vector2);
    if (vector1.length !== vector2.length)
        throw 'Vectors do not have the same dimensions';
    let cov = 0;
    const l = vector1.length;
    for (let i = 0; i < l; i++) {
        let x = vector1[i] - mean1;
        let y = vector2[i] - mean2;
        cov += x * y;
    }
    if (unbiased)
        return cov / (l - 1);
    else
        return cov / l;
}
/**
 *
 * @param {number[]} values
 * @param unbiased = true
 * @returns
 */
export function skewness(values, unbiased = true) {
    let theMean = mean(values);
    let s2 = 0, s3 = 0;
    const l = values.length;
    for (let i = 0; i < l; i++) {
        let dev = values[i] - theMean;
        s2 += dev * dev;
        s3 += dev * dev * dev;
    }
    let m2 = s2 / l;
    let m3 = s3 / l;
    let g = m3 / (Math.pow(m2, 3 / 2.0));
    if (unbiased) {
        let a = Math.sqrt(l * (l - 1));
        let b = l - 2;
        return (a / b) * g;
    }
    else {
        return g;
    }
}
/**
 *
 * @param {number[]} values
 * @param unbiased = true
 * @returns
 */
export function kurtosis(values, unbiased = true) {
    let theMean = mean(values);
    const n = values.length;
    let s2 = 0, s4 = 0;
    for (let i = 0; i < n; i++) {
        let dev = values[i] - theMean;
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
        return a * b - 3 * c;
    }
    else {
        return m4 / (m2 * m2) - 3;
    }
}
/**
 *
 * @param {number[]} values
 * @param eps = 0
 * @returns
 */
export function entropy(values, eps = 0) {
    let sum = 0;
    let l = values.length;
    for (let i = 0; i < l; i++)
        sum += values[i] * Math.log(values[i] + eps);
    return -sum;
}
/**
 *
 * @param {number[]} values
 * @param {number[]} weights
 * @returns
 */
export function weightedMean(values, weights) {
    let sum = 0;
    let l = values.length;
    for (let i = 0; i < l; i++)
        sum += values[i] * weights[i];
    return sum;
}
/**
 *
 * @param {number[]} values
 * @param {number[]} weights
 * @returns
 */
export function weightedStandardDeviation(values, weights) {
    return Math.sqrt(weightedVariance(values, weights));
}
/**
 *
 * @param {number[]} values
 * @param {number[]} weights
 * @returns
 */
export function weightedVariance(values, weights) {
    let theMean = weightedMean(values, weights);
    let vari = 0;
    let l = values.length;
    let a = 0, b = 0;
    for (let i = 0; i < l; i++) {
        let z = values[i] - theMean;
        let w = weights[i];
        vari += w * (z * z);
        b += w;
        a += w * w;
    }
    return vari * (b / (b * b - a));
}
/**
 *
 * @param {number[]} values
 * @param inPlace = false
 * @returns
 */
export function center(values, inPlace = false) {
    let result = values;
    if (!inPlace)
        result = [].concat(values);
    let theMean = mean(result);
    let l = result.length;
    for (let i = 0; i < l; i++) {
        result[i] -= theMean;
    }
    return result; //added return
}
/**
 *
 * @param {number[]} values
 * @param standardDev = standardDeviation(values)
 * @param inPlace = false
 * @returns
 */
export function standardize(values, standardDev = standardDeviation(values), inPlace = false) {
    let l = values.length;
    let result = inPlace ? values : new Array(l);
    for (let i = 0; i < l; i++)
        result[i] = values[i] / standardDev;
    return result;
}
/**
 *
 * @param {number[]} array
 * @returns
 */
export function cumulativeSum(array) {
    let l = array.length;
    let result = new Array(l);
    result[0] = array[0];
    for (let i = 1; i < l; i++)
        result[i] = result[i - 1] + array[i];
    return result;
}
//# sourceMappingURL=array.js.map