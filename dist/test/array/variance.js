'use strict';
import { array as stat } from '../../src/index';
describe('variance and standard deviation', () => {
    let data = [15, 13, 17, 7];
    it('variance', () => {
        let v = stat.variance(data);
        v.should.be.approximately(18.667, 1e-3);
        stat.variance(data, true).should.equal(v);
        stat.variance(data, false).should.equal(14);
    });
    it('standard deviation', () => {
        let s = stat.standardDeviation(data);
        stat.standardDeviation(data, true).should.equal(s);
        stat.standardDeviation(data, false).should.equal(Math.sqrt(14));
    });
});
//# sourceMappingURL=variance.js.map