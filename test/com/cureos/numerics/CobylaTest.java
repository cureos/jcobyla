/*
 * jcobyla
 * 
 * The MIT License
 *
 * Copyright (c) 2012 Anders Gustafsson, Cureos AB.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files 
 * (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, 
 * publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, 
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
 * FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION 
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 * 
 * Remarks:
 * 
 * The original Fortran 77 version of this code was by Michael Powell (M.J.D.Powell @ damtp.cam.ac.uk)
 * The Fortran 90 version was by Alan Miller (Alan.Miller @ vic.cmis.csiro.au). Latest revision - 30 October 1998
 */
package com.cureos.numerics;

import static org.junit.Assert.assertArrayEquals;
import org.junit.Test;

/**
 * Test class for COBYLA2 employing tests from Report DAMTP 1992/NA5.
 * 
 * @author Anders Gustafsson, Cureos AB.
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
        System.out.format("%nOutput from test problem 1 (Simple quadratic)%n");
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
        System.out.format("%nOutput from test problem 2 (2D unit circle calculation)%n");
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
        System.out.format("%nOutput from test problem 3 (3D ellipsoid calculation)%n");
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
     * Weak version of Rosenbrock's problem.
     */
    @Test
    public void test04FindMinimum() {
        System.out.format("%nOutput from test problem 4 (Weak Rosenbrock)%n");
        Calcfc calcfc = new Calcfc() {
            @Override
            public double Compute(int n, int m, double[] x, double[] con) {
                return Math.pow(x[0] * x[0] - x[1], 2.0) + Math.pow(1.0 + x[0], 2.0);
            }
        };
        double[] x = {1.0, 1.0 };
        CobylaExitStatus result = Cobyla.FindMinimum(calcfc, 2, 0, x, rhobeg, rhoend, iprint, maxfun);
        assertArrayEquals(null, new double[] { -1.0, 1.0 }, x, 1.0e-4);
    }
    
    /**
     * Intermediate version of Rosenbrock's problem.
     */
    @Test
    public void test05FindMinimum() {
        System.out.format("%nOutput from test problem 5 (Intermediate Rosenbrock)%n");
        Calcfc calcfc = new Calcfc() {
            @Override
            public double Compute(int n, int m, double[] x, double[] con) {
                return 10.0 * Math.pow(x[0] * x[0] - x[1], 2.0) + Math.pow(1.0 + x[0], 2.0);
            }
        };
        double[] x = {1.0, 1.0 };
        CobylaExitStatus result = Cobyla.FindMinimum(calcfc, 2, 0, x, rhobeg, rhoend, iprint, maxfun);
        assertArrayEquals(null, new double[] { -1.0, 1.0 }, x, 3.0e-4);
    }
    
    /**
     * This problem is taken from Fletcher's book Practical Methods of
     * Optimization and has the equation number (9.1.15).
     */
    @Test
    public void test06FindMinimum() {
        System.out.format("%nOutput from test problem 6 (Equation (9.1.15) in Fletcher's book)%n");
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
    
    /**
     * This problem is taken from Fletcher's book Practical Methods of
     * Optimization and has the equation number (14.4.2).
     */
    @Test
    public void test07FindMinimum() {
        System.out.format("%nOutput from test problem 7 (Equation (14.4.2) in Fletcher)%n");
        Calcfc calcfc = new Calcfc() {
            @Override
            public double Compute(int n, int m, double[] x, double[] con) {
                con[0] = 5.0 * x[0] - x[1] + x[2];
                con[1] = x[2] - x[0] * x[0] - x[1] * x[1] - 4.0 * x[1];
                con[2] = x[2] - 5.0 * x[0] - x[1];
                return x[2];
            }
        };
        double[] x = {1.0, 1.0, 1.0 };
        CobylaExitStatus result = Cobyla.FindMinimum(calcfc, 3, 3, x, rhobeg, rhoend, iprint, maxfun);
        assertArrayEquals(null, new double[] { 0.0, -3.0, -3.0 }, x, 1.0e-5);
    }
    
    /**
     * This problem is taken from page 66 of Hock and Schittkowski's book Test
     * Examples for Nonlinear Programming Codes. It is their test problem Number
     * 43, and has the name Rosen-Suzuki.
     */
    @Test
    public void test08FindMinimum() {
        System.out.format("%nOutput from test problem 8 (Rosen-Suzuki)%n");
        Calcfc calcfc = new Calcfc() {
            @Override
            public double Compute(int n, int m, double[] x, double[] con) {
                con[0] = 8.0 - x[0] * x[0] - x[1] * x[1] - x[2] * x[2] - x[3] * x[3] - x[0] + x[1] - x[2] + x[3];
                con[1] = 10.0 - x[0] * x[0] - 2.0 * x[1] * x[1] - x[2] * x[2] - 2.0 * x[3] * x[3] + x[0] + x[3];
                con[2] = 5.0 - 2.0 * x[0] * x[0] - x[1] * x[1] - x[2] * x[2] - 2.0 * x[0] + x[1] + x[3];
                return x[0] * x[0] + x[1] * x[1] + 2.0 * x[2] * x[2] + x[3] * x[3] - 5.0 * x[0] - 
                        5.0 * x[1] - 21.0 * x[2] + 7.0 * x[3];
            }
        };
        double[] x = {1.0, 1.0, 1.0, 1.0 };
        CobylaExitStatus result = Cobyla.FindMinimum(calcfc, 4, 3, x, rhobeg, rhoend, iprint, maxfun);
        assertArrayEquals(null, new double[] { 0.0, 1.0, 2.0, -1.0 }, x, 1.0e-5);
    }
    
    /**
     * This problem is taken from page 111 of Hock and Schittkowski's
     * book Test Examples for Nonlinear Programming Codes. It is their
     * test problem Number 100.
     */
    @Test
    public void test09FindMinimum() {
        System.out.format("%nOutput from test problem 9 (Hock and Schittkowski 100)%n");
        Calcfc calcfc = new Calcfc() {
            @Override
            public double Compute(int n, int m, double[] x, double[] con) {
                con[0] = 127.0 - 2.0 * x[0] * x[0] - 3.0 * Math.pow(x[1], 4.0) - x[2] - 4.0 * x[3] * x[3] - 5.0 * x[4];
                con[1] = 282.0 - 7.0 * x[0] - 3.0 * x[1] - 10.0 * x[2] * x[2] - x[3] + x[4];
                con[2] = 196.0 - 23.0 * x[0] - x[1] * x[1] - 6.0 * x[5] * x[5] + 8.0 * x[6];
                con[3] = -4.0 * x[0] * x[0] - x[1] * x[1] + 3.0 * x[0] * x[1] - 2.0 * x[2] * x[2] - 5.0 * x[5] + 11.0 * x[6];
                return Math.pow(x[0] - 10.0, 2.0) + 5.0 * Math.pow(x[1] - 12.0, 2.0) + Math.pow(x[2], 4.0) +
                        3.0 * Math.pow(x[3] - 11.0, 2.0) + 10.0 * Math.pow(x[4], 6.0) + 7.0 * x[5] * x[5] + Math.pow(x[6], 4.0) -
                        4.0 * x[5] * x[6] - 10.0 * x[5] - 8.0 * x[6];
            }
        };
        double[] x = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
        CobylaExitStatus result = Cobyla.FindMinimum(calcfc, 7, 4, x, rhobeg, rhoend, iprint, maxfun);
        assertArrayEquals(null, 
                new double[] { 2.330499, 1.951372, -0.4775414, 4.365726, -0.624487, 1.038131, 1.594227 }, x, 1.0e-5);
    }
    
    /**
     * This problem is taken from page 415 of Luenberger's book Applied 
     * Nonlinear Programming. It is to maximize the area of a hexagon of
     * unit diameter.
     */
    @Test
    public void test10FindMinimum() {
        System.out.format("%nOutput from test problem 10 (Hexagon area)%n");
        Calcfc calcfc = new Calcfc() {
            @Override
            public double Compute(int n, int m, double[] x, double[] con) {
                con[0] = 1.0 - x[2] * x[2] - x[3] * x[3];
                con[1] = 1.0 - x[8] * x[8];
                con[2] = 1.0 - x[4] * x[4] - x[5] * x[5];
                con[3] = 1.0 - x[0] * x[0] - Math.pow(x[1] - x[8], 2.0);
                con[4] = 1.0 - Math.pow(x[0] - x[4], 2.0) - Math.pow(x[1] - x[5], 2.0);
                con[5] = 1.0 - Math.pow(x[0] - x[6], 2.0) - Math.pow(x[1] - x[7], 2.0);
                con[6] = 1.0 - Math.pow(x[2] - x[4], 2.0) - Math.pow(x[3] - x[5], 2.0);
                con[7] = 1.0 - Math.pow(x[2] - x[6], 2.0) - Math.pow(x[3] - x[7], 2.0);
                con[8] = 1.0 - x[6] * x[6] - Math.pow(x[7] - x[8], 2.0);
                con[9] = x[0] * x[3] - x[1] * x[2];
                con[10] = x[2] * x[8];
                con[11] = -x[4] * x[8];
                con[12] = x[4] * x[7] - x[5] * x[6];
                con[13] = x[8];
                return -0.5 * (x[0] * x[3] - x[1] * x[2] + x[2] * x[8] - x[4] * x[8] + x[4] * x[7] - x[5] * x[6]);
            }
        };
        double[] x = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
        CobylaExitStatus result = Cobyla.FindMinimum(calcfc, 9, 14, x, rhobeg, rhoend, iprint, maxfun);
        assertArrayEquals(null, 
                new double[] { x[0], x[1], x[2], x[3], x[0], x[1], x[2], x[3], 0.0 }, x, 1.0e-4);
    }
}
