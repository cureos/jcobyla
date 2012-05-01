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
    
    // FIELDS
    
    private double rhobeg = 0.5;
    private double rhoend = 1.0e-6;
    private int iprint = 1;
    private int maxfun = 3500;

    // TESTS
    
    /**
     * Minimization of a simple quadratic function of two variables.
     */
    @Test
    public void test01FindMinimum() {
        System.out.format("%n1 (Simple quadratic)%n");
        Calcfc calcfc = new Calcfc() {
            @Override
            public double Compute(int n, int m, double[] x, double[] con) {
                return 10.0 * Math.pow(x[0] + 1.0, 2.0) + Math.pow(x[1], 2.0);
            }
        };
        double[] x = {1.0, 1.0 };
        CobylaExitStatus result = Cobyla.FindMinimum(calcfc, 2, 0, x, rhobeg, rhoend, iprint, maxfun);
        assertArrayEquals(null, new double[] { -1.0, 0.0 }, x, 1.0e-5);
    }
    
    /**
     * Easy two dimensional minimization in unit circle.
     */
    @Test
    public void test02FindMinimum() {
        System.out.format("%n2 (2D unit circle calculation)%n");
        Calcfc calcfc = new Calcfc() {
            @Override
            public double Compute(int n, int m, double[] x, double[] con) {
                con[0] = 1.0 - x[0] * x[0] - x[1] * x[1];
                return x[0] * x[1];
            }
        };
        double[] x = {1.0, 1.0 };
        CobylaExitStatus result = Cobyla.FindMinimum(calcfc, 2, 1, x, rhobeg, rhoend, iprint, maxfun);
        assertArrayEquals(null, new double[] { Math.sqrt(0.5), -Math.sqrt(0.5) }, x, 1.0e-5);
    }
    
    /**
     * Easy three dimensional minimization in ellipsoid.
     */
    @Test
    public void test03FindMinimum() {
        System.out.format("%n3 (3D ellipsoid calculation)%n");
        Calcfc calcfc = new Calcfc() {
            @Override
            public double Compute(int n, int m, double[] x, double[] con) {
                con[0] = 1.0 - x[0] * x[0] - 2.0 * x[1] * x[1] - 3.0 * x[2] * x[2];
                return x[0] * x[1] * x[2];
            }
        };
        double[] x = {1.0, 1.0, 1.0 };
        CobylaExitStatus result = Cobyla.FindMinimum(calcfc, 3, 1, x, rhobeg, rhoend, iprint, maxfun);
        assertArrayEquals(null, new double[] { 1.0 / Math.sqrt(3.0), 1.0 / Math.sqrt(6.0), -1.0 / 3.0 }, x, 1.0e-5);
    }
    
    /**
     * This problem is taken from Fletcher's book Practical Methods of
     * Optimization and has the equation number (9.1.15).
     */
    @Test
    public void test06FindMinimum() {
        System.out.format("%n6 (Equation (9.1.15) in Fletcher's book)%n");
        Calcfc calcfc = new Calcfc() {
            @Override
            public double Compute(int n, int m, double[] x, double[] con) {
                con[0] = x[1] - x[0] * x[0];
                con[1] = 1.0 - x[0] * x[0] - x[1] * x[1];
                return -x[0] - x[1];
            }
        };
        double[] x = {1.0, 1.0 };
        CobylaExitStatus result = Cobyla.FindMinimum(calcfc, 2, 2, x, rhobeg, rhoend, iprint, maxfun);
        assertArrayEquals(null, new double[] { Math.sqrt(0.5), Math.sqrt(0.5) }, x, 1.0e-5);
    }
}
