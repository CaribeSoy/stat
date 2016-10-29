'use strict';
import { array as stat } from '../../src/index';
describe('min-max', () => {
    let arr = [5, 10, -4, 3, 1, -16, 32, 14, -3, 9];
    it('should find the min value', () => {
        stat.min(arr).should.equal(-16);
    });
    it('should find the max value', () => {
        stat.max(arr).should.equal(32);
    });
    it('should find min and max values', () => {
        let result = stat.minMax(arr);
        result.should.be.an.Object();
        result.min.should.equal(-16);
        result.max.should.equal(32);
    });
});
//# sourceMappingURL=min-max.js.map