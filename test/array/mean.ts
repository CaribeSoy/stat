'use strict';

import * as stat from '../../src/array';
const should = require('should');

describe('mean', () => {

    let arr1 = [4, 36, 45, 50, 75];
    let arr2 = [3, 3, 3];
    let unsorted = [3, 1, 6, 14, 12];
    let sorted = [1, 3, 6, 12, 14];

    it('sum', () => {
        stat.sum(arr1).should.equal(210);
        stat.sum(arr2).should.equal(9);
    });

    it('arithmetic mean', () => {
        stat.mean(arr1).should.equal(42);
        stat.arithmeticMean(arr1).should.equal(42);
        stat.mean(arr2).should.equal(3);
    });

    it('geometric mean', () => {
        stat.geometricMean(arr1).should.be.approximately(30, 1e-6);
        stat.geometricMean(arr2).should.equal(3);
    });

    it('log mean', () => {
        Math.exp(stat.logMean(arr1)).should.be.approximately(30, 1e-6);
    });

    it('grand mean', () => {
        let means = [1, 2, 3];
        let sizes = [2, 5, 3];
        stat.grandMean(means, sizes).should.equal(2.1);
    });

    it('truncated mean', () => {
        stat.truncatedMean(unsorted, 0.2).should.equal(7);
        stat.truncatedMean(sorted, 0.2, true).should.equal(7);
        stat.truncatedMean(sorted, 0.5).should.equal(6);
    });

    it('harmonic mean', () => {
        let arr1 = [1, 2, 4];
        let arr2 = [3, 0, 6];
        stat.harmonicMean(arr1).should.be.approximately(1.714, 1e-3);
        (() => {
            stat.harmonicMean(arr2);
        }).should.throw(RangeError);
    });

    it('contraharmonic mean', () => {
        let arr1 = [1, 5, 6];
        let arr2 = [2, -6, -1];
        stat.contraHarmonicMean(arr1).should.be.approximately(5.167, 1e-3);
        (() => {
            stat.contraHarmonicMean(arr2);
        }).should.throw(RangeError);
    });

    it('median', () => {
        stat.median(unsorted).should.equal(6);
        stat.median(sorted, true).should.equal(6);
        stat.median([2, 4, 6, 8]).should.equal(5);
    });

});