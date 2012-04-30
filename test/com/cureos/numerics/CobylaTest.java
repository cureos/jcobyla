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
     * Test 2 of FindMinimum method, of class Cobyla.
     */
    @Test
    public void test2FindMinimum() {
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
        assertArrayEquals("2 (2D unit circle calculation)", 
                new double[] { Math.sqrt(0.5), -Math.sqrt(0.5) }, x, 1.0e-5);
    }
}
