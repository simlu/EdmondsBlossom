import org.junit.Assert;
import org.junit.Test;


public class BlossomTest {
    @Test
    public void test10_empty() {
        // empty input graph
        Assert.assertArrayEquals(
                new Blossom(new int[0][], false).maxWeightMatching(),
                new int[0]
        );
    }

    @Test
    public void test11_singleedge() {
        // single edge
        Assert.assertArrayEquals(
                new Blossom(new int[][]{new int[]{0, 1, 1}}, false).maxWeightMatching(),
                new int[]{1, 0}
        );
    }

    @Test
    public void test12() {
        Assert.assertArrayEquals(
                new Blossom(new int[][]{new int[]{1, 2, 10}, new int[]{2, 3, 11}}, false).maxWeightMatching(),
                new int[]{-1, -1, 3, 2}
        );
    }

    @Test
    public void test13() {
        Assert.assertArrayEquals(
                new Blossom(new int[][]{new int[]{1, 2, 5}, new int[]{2, 3, 11}, new int[]{3, 4, 5}}, false).maxWeightMatching(),
                new int[]{-1, -1, 3, 2, -1}
        );
    }

    @Test
    public void test14_maxcard() {
        // maximum cardinality
        Assert.assertArrayEquals(
                new Blossom(new int[][]{new int[]{1, 2, 5}, new int[]{2, 3, 11}, new int[]{3, 4, 5}}, true).maxWeightMatching(),
                new int[]{-1, 2, 1, 4, 3}
        );
    }

    @Test
    public void test15_float() {
        // floating point weigths
        Assert.assertArrayEquals(
                new Blossom(new float[][]{new float[]{1, 2, (float) Math.PI}, new float[]{2, 3, (float) Math.exp(1)}, new float[]{1, 3, 3.0f}, new float[]{1, 4, (float) Math.sqrt(2.0)}}, false).maxWeightMatching(),
                new int[]{-1, 4, 3, 2, 1}
        );
    }

    @Test
    public void test16_negative() {
        // negative weights
        Assert.assertArrayEquals(
                new Blossom(new int[][]{new int[]{1, 2, 2}, new int[]{1, 3, -2}, new int[]{2, 3, 1}, new int[]{2, 4, -1}, new int[]{3, 4, -6}}, false).maxWeightMatching(),
                new int[]{-1, 2, 1, -1, -1}
        );
        Assert.assertArrayEquals(
                new Blossom(new int[][]{new int[]{1, 2, 2}, new int[]{1, 3, -2}, new int[]{2, 3, 1}, new int[]{2, 4, -1}, new int[]{3, 4, -6}}, true).maxWeightMatching(),
                new int[]{-1, 3, 4, 1, 2}
        );
    }

    @Test
    public void test20_sblossom() {
        // create S-blossom and use it for augmentation
        Assert.assertArrayEquals(
                new Blossom(new int[][]{new int[]{1, 2, 8}, new int[]{1, 3, 9}, new int[]{2, 3, 10}, new int[]{3, 4, 7}}, false).maxWeightMatching(),
                new int[]{-1, 2, 1, 4, 3}
        );
        Assert.assertArrayEquals(
                new Blossom(new int[][]{new int[]{1, 2, 8}, new int[]{1, 3, 9}, new int[]{2, 3, 10}, new int[]{3, 4, 7}, new int[]{1, 6, 5}, new int[]{4, 5, 6}}, false).maxWeightMatching(),
                new int[]{-1, 6, 3, 2, 5, 4, 1}
        );
    }

    @Test
    public void test21_tblossom() {
        // create S-blossom, relabel as T-blossom, use for augmentation
        Assert.assertArrayEquals(
                new Blossom(new int[][]{new int[]{1, 2, 9}, new int[]{1, 3, 8}, new int[]{2, 3, 10}, new int[]{1, 4, 5}, new int[]{4, 5, 4}, new int[]{1, 6, 3}}, false).maxWeightMatching(),
                new int[]{-1, 6, 3, 2, 5, 4, 1}
        );
        Assert.assertArrayEquals(
                new Blossom(new int[][]{new int[]{1, 2, 9}, new int[]{1, 3, 8}, new int[]{2, 3, 10}, new int[]{1, 4, 5}, new int[]{4, 5, 3}, new int[]{1, 6, 4}}, false).maxWeightMatching(),
                new int[]{-1, 6, 3, 2, 5, 4, 1}
        );
        Assert.assertArrayEquals(
                new Blossom(new int[][]{new int[]{1, 2, 9}, new int[]{1, 3, 8}, new int[]{2, 3, 10}, new int[]{1, 4, 5}, new int[]{4, 5, 3}, new int[]{3, 6, 4}}, false).maxWeightMatching(),
                new int[]{-1, 2, 1, 6, 5, 4, 3}
        );
    }

    @Test
    public void test22_s_nest() {
        // create nested S-blossom, use for augmentation
        Assert.assertArrayEquals(
                new Blossom(new int[][]{new int[]{1, 2, 9}, new int[]{1, 3, 9}, new int[]{2, 3, 10}, new int[]{2, 4, 8}, new int[]{3, 5, 8}, new int[]{4, 5, 10}, new int[]{5, 6, 6}}, false).maxWeightMatching(),
                new int[]{-1, 3, 4, 1, 2, 6, 5}
        );
    }

    @Test
    public void test23_s_relabel_nest() {
        // create S-blossom, relabel as S, include in nested S-blossom
        Assert.assertArrayEquals(
                new Blossom(new int[][]{new int[]{1, 2, 10}, new int[]{1, 7, 10}, new int[]{2, 3, 12}, new int[]{3, 4, 20}, new int[]{3, 5, 20}, new int[]{4, 5, 25}, new int[]{5, 6, 10}, new int[]{6, 7, 10}, new int[]{7, 8, 8}}, false).maxWeightMatching(),
                new int[]{-1, 2, 1, 4, 3, 6, 5, 8, 7}
        );
    }

    @Test
    public void test24_s_nest_expand() {
        // create nested S-blossom, augment, expand recursively
        Assert.assertArrayEquals(
                new Blossom(new int[][]{new int[]{1, 2, 8}, new int[]{1, 3, 8}, new int[]{2, 3, 10}, new int[]{2, 4, 12}, new int[]{3, 5, 12}, new int[]{4, 5, 14}, new int[]{4, 6, 12}, new int[]{5, 7, 12}, new int[]{6, 7, 14}, new int[]{7, 8, 12}}, false).maxWeightMatching(),
                new int[]{-1, 2, 1, 5, 6, 3, 4, 8, 7}
        );
    }

    @Test
    public void test25_s_t_expand() {
        // create S-blossom, relabel as T, expand
        Assert.assertArrayEquals(
                new Blossom(new int[][]{new int[]{1, 2, 23}, new int[]{1, 5, 22}, new int[]{1, 6, 15}, new int[]{2, 3, 25}, new int[]{3, 4, 22}, new int[]{4, 5, 25}, new int[]{4, 8, 14}, new int[]{5, 7, 13}}, false).maxWeightMatching(),
                new int[]{-1, 6, 3, 2, 8, 7, 1, 5, 4}
        );
    }

    @Test
    public void test26_s_nest_t_expand() {
        // create nested S-blossom, relabel as T, expand
        Assert.assertArrayEquals(
                new Blossom(new int[][]{new int[]{1, 2, 19}, new int[]{1, 3, 20}, new int[]{1, 8, 8}, new int[]{2, 3, 25}, new int[]{2, 4, 18}, new int[]{3, 5, 18}, new int[]{4, 5, 13}, new int[]{4, 7, 7}, new int[]{5, 6, 7}}, false).maxWeightMatching(),
                new int[]{-1, 8, 3, 2, 7, 6, 5, 4, 1}
        );
    }

    @Test
    public void test30_tnasty_expand() {
        // create blossom, relabel as T in more than one way, expand, augment
        Assert.assertArrayEquals(
                new Blossom(new int[][]{new int[]{1, 2, 45}, new int[]{1, 5, 45}, new int[]{2, 3, 50}, new int[]{3, 4, 45}, new int[]{4, 5, 50}, new int[]{1, 6, 30}, new int[]{3, 9, 35}, new int[]{4, 8, 35}, new int[]{5, 7, 26}, new int[]{9, 10, 5}}, false).maxWeightMatching(),
                new int[]{-1, 6, 3, 2, 8, 7, 1, 5, 4, 10, 9}
        );
    }

    @Test
    public void test31_tnasty2_expand() {
        // again but slightly different
        Assert.assertArrayEquals(
                new Blossom(new int[][]{new int[]{1, 2, 45}, new int[]{1, 5, 45}, new int[]{2, 3, 50}, new int[]{3, 4, 45}, new int[]{4, 5, 50}, new int[]{1, 6, 30}, new int[]{3, 9, 35}, new int[]{4, 8, 26}, new int[]{5, 7, 40}, new int[]{9, 10, 5}}, false).maxWeightMatching(),
                new int[]{-1, 6, 3, 2, 8, 7, 1, 5, 4, 10, 9}
        );
    }

    @Test
    public void test32_t_expand_leastslack() {
        // create blossom, relabel as T, expand such that a new least-slack S-to-free edge is produced, augment
        Assert.assertArrayEquals(
                new Blossom(new int[][]{new int[]{1, 2, 45}, new int[]{1, 5, 45}, new int[]{2, 3, 50}, new int[]{3, 4, 45}, new int[]{4, 5, 50}, new int[]{1, 6, 30}, new int[]{3, 9, 35}, new int[]{4, 8, 28}, new int[]{5, 7, 26}, new int[]{9, 10, 5}}, false).maxWeightMatching(),
                new int[]{-1, 6, 3, 2, 8, 7, 1, 5, 4, 10, 9}
        );
    }

    @Test
    public void test33_nest_tnasty_expand() {
        // create nested blossom, relabel as T in more than one way, expand outer blossom such that inner blossom ends up on an augmenting path
        Assert.assertArrayEquals(
                new Blossom(new int[][]{new int[]{1, 2, 45}, new int[]{1, 7, 45}, new int[]{2, 3, 50}, new int[]{3, 4, 45}, new int[]{4, 5, 95}, new int[]{4, 6, 94}, new int[]{5, 6, 94}, new int[]{6, 7, 50}, new int[]{1, 8, 30}, new int[]{3, 11, 35}, new int[]{5, 9, 36}, new int[]{7, 10, 26}, new int[]{11, 12, 5}}, false).maxWeightMatching(),
                new int[]{-1, 8, 3, 2, 6, 9, 4, 10, 1, 5, 7, 12, 11}
        );
    }

    @Test
    public void test34_nest_relabel_expand() {
        // create nested S-blossom, relabel as S, expand recursively
        Assert.assertArrayEquals(
                new Blossom(new int[][]{new int[]{1, 2, 40}, new int[]{1, 3, 40}, new int[]{2, 3, 60}, new int[]{2, 4, 55}, new int[]{3, 5, 55}, new int[]{4, 5, 50}, new int[]{1, 8, 15}, new int[]{5, 7, 30}, new int[]{7, 6, 10}, new int[]{8, 10, 10}, new int[]{4, 9, 30}}, false).maxWeightMatching(),
                new int[]{-1, 2, 1, 5, 9, 3, 7, 6, 10, 4, 8}
        );
    }

    // below is ported from: https://github.com/mattkrick/EdmondsBlossom/blob/master/spec/blossomSpec.js

    @Test
    public void test_other_int() {
        Assert.assertArrayEquals(
                new Blossom(new int[][]{new int[]{1, 2, 5}, new int[]{1, 3, 6}, new int[]{1, 4, 2}, new int[]{2, 3, 4}, new int[]{2, 4, 1}, new int[]{3, 4, 3}}, false).maxWeightMatching(),
                new int[]{-1, 2, 1, 4, 3}
        );
    }

    @Test
    public void test_other_float() {
        Assert.assertArrayEquals(
                new Blossom(new float[][]{new float[]{1, 2, 47.2612f}, new float[]{1, 3, 46.9176f}, new float[]{2, 3, 49.3305f}, new float[]{1, 4, 44.7978f}, new float[]{2, 4, 49.1123f}, new float[]{2, 5, 51.1539f}, new float[]{4, 5, 50.5430f}, new float[]{2, 6, 48.2873f}, new float[]{3, 6, 47.7470f}, new float[]{4, 6, 46.8674f}, new float[]{5, 6, 48.8397f}}, false).maxWeightMatching(),
                new int[]{-1, 3, 6, 1, 5, 4, 2}
        );
    }

    @Test
    public void test_simple_zero_index_love_triangle() {
        // should handle a simple zero-index love triangle
        Assert.assertArrayEquals(
                new Blossom(new int[][]{new int[]{0, 1, 6}, new int[]{0, 2, 10}, new int[]{1, 2, 5}}, false).maxWeightMatching(),
                new int[]{2, -1, 0}
        );
    }

    @Test
    public void test_online_example1() {
        // Reference: http://uoj.ac/problem/81
        Assert.assertArrayEquals(
                new Blossom(new int[][] {new int[]{4, 6, 9}, new int[]{2, 6, 4}, new int[]{2, 5, 6}, new int[]{1, 4, 8}, new int[]{4, 0, 9}, new int[]{0, 2, 6}, new int[]{5, 4, 1}, new int[]{1, 6, 4}, new int[]{1, 2, 5}, new int[]{5, 3, 2}, new int[]{6, 0, 5}, new int[]{4, 3, 4}, new int[]{3, 0, 3}, new int[]{4, 2, 9}, new int[]{6, 5, 4}, new int[]{1, 0, 3}, new int[]{3, 2, 9}, new int[]{5, 1, 7}, new int[]{3, 1, 8}, new int[]{5, 0, 10}}, true).maxWeightMatching(),
                new int[]{5, -1, 3, 2, 6, 0, 4}
        );
    }

}