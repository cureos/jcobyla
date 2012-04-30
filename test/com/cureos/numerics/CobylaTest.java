/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cureos.numerics;

import org.junit.*;
import static org.junit.Assert.*;

/**
 *
 * @author anders
 */
public class CobylaTest {
    
    /**
     * Easy two dimensional minimization in unit circle.
     */
    @Test
    public void test02FindMinimum() {
        System.out.println("2 (2D unit circle calculation)");
        Calcfc calcfc = new Calcfc() {
            @Override
            public double Compute(int n, int m, double[] x, double[] con) {
                con[0] = 1.0 - x[0] * x[0] - x[1] * x[1];
                return x[0] * x[1];
            }
        };
        double[] x = {1.0, 1.0 };
        CobylaExitStatus result = Cobyla.FindMinimum(calcfc, 2, 1, x, 0.5, 1.0e-6, 2, 3500);
        assertArrayEquals(null, new double[] { Math.sqrt(0.5), -Math.sqrt(0.5) }, x, 1.0e-5);
    }
    
    /**
     * This problem is taken from Fletcher's book Practical Methods of
     * Optimization and has the equation number (9.1.15).
     */
    @Test
    public void test06FindMinimum() {
        System.out.println("6 (Equation (9.1.15) in Fletcher's book)");
        Calcfc calcfc = new Calcfc() {
            @Override
            public double Compute(int n, int m, double[] x, double[] con) {
                con[0] = x[1] - x[0] * x[0];
                con[1] = 1.0 - x[0] * x[0] - x[1] * x[1];
                return -x[0] - x[1];
            }
        };
        double[] x = {1.0, 1.0 };
        CobylaExitStatus result = Cobyla.FindMinimum(calcfc, 2, 2, x, 0.5, 1.0e-6, 2, 3500);
        assertArrayEquals(null, new double[] { Math.sqrt(0.5), Math.sqrt(0.5) }, x, 1.0e-5);
    }
}
