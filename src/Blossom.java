/* Converted to Java from Python by Matt Krick. Original: http://jorisvr.nl/maximummatching.html */
// Javascript reference: https://github.com/mattkrick/EdmondsBlossom

import it.unimi.dsi.fastutil.ints.IntArrayList;

import java.util.Arrays;
import java.util.Collections;

public class Blossom {
    private static final boolean CHECK_DELTA = false;
    private static final boolean CHECK_OPTIMUM = false;
    private static final boolean DEBUG = false;
    private static final Tracker tracker = null;

    private final int[][] edges;
    private final boolean maxcardinality;
    private final int nedge;
    private final int nvertex;
    private final int maxweight;
    private final int[] endpoint;
    private final IntArrayList[] neighbend;
    private final int[] mate;
    private final int[] label;
    private final int[] labelend;
    private final int[] inblossom;
    private final int[] blossomparent;
    private final IntArrayList[] blossomchilds;
    private final int[] blossombase;
    private final IntArrayList[] blossomendps;
    private final int[] bestedge;
    private final int[][] blossombestedges;
    private final IntArrayList unusedblossoms;
    private final int[] dualvar;
    private final boolean[] allowedge;
    private final IntArrayList queue;

    private static int[][] getIntEdges(float[][] edges) {
        // Convert float edge weights to integers
        int[][] result = new int[edges.length][];
        // find best scaling value
        int safeScaleValue = Integer.MAX_VALUE;
        for (float[] edge : edges) {
            safeScaleValue = Math.min(safeScaleValue, (int) Math.floor(Integer.MAX_VALUE / Math.abs(edge[2])));
        }
        safeScaleValue = (int) Math.pow(10, Integer.toString(safeScaleValue).length() - 1d);

        // convert floats to integers
        for (int i = 0; i < edges.length; i++) {
            float[] edge = edges[i];
            assert Float.floatToIntBits(edge[0]) == Float.floatToIntBits(Math.round(edge[0]));
            assert Float.floatToIntBits(edge[1]) == Float.floatToIntBits(Math.round(edge[1]));
            result[i] = new int[]{(int) edge[0], (int) edge[1], (int) (safeScaleValue * edge[2])};
        }
        return result;
    }

    public Blossom(float[][] edges, boolean maxcardinality) {
        // Helper to allow for float weights
        this(getIntEdges(edges), maxcardinality);
    }

    public Blossom(int[][] edges, boolean maxcardinality) {
        // Vertices are numbered 0 .. (nvertex-1).
        // Non-trivial blossoms are numbered nvertex .. (2*nvertex-1)
        //
        // Edges are numbered 0 .. (nedge-1).
        // Edge endpoints are numbered 0 .. (2*nedge-1), such that endpoints
        // (2*k) and (2*k+1) both belong to edge k.
        //
        // Many terms used in the comments (sub-blossom, T-vertex) come from
        // the paper by Galil; read the paper before reading this code.
        this.edges = edges;
        this.maxcardinality = maxcardinality;

        // Count vertices.
        int nedge = edges.length;
        int nvertex = 0;
        for (int[] edge : edges) {
            assert edge[0] >= 0;
            assert edge[1] >= 0;
            assert edge[0] != edge[1];
            if (edge[0] >= nvertex) {
                nvertex = edge[0] + 1;
            }
            if (edge[1] >= nvertex) {
                nvertex = edge[1] + 1;
            }
        }
        this.nedge = nedge;
        this.nvertex = nvertex;

        // Find the maximum edge weight.
        int maxweight = 0;
        for (int[] edge : edges) {
            maxweight = Math.max(edge[2], maxweight);
        }
        this.maxweight = maxweight;

        // If p is an edge endpoint,
        // endpoint[p] is the vertex to which endpoint p is attached.
        // Not modified by the algorithm.
        int[] endpoint = new int[2 * nedge];
        for (int p = 0; p < 2 * nedge; p++) {
            endpoint[p] = edges[p / 2][p % 2];
        }
        this.endpoint = endpoint;

        // If v is a vertex,
        // neighbend[v] is the list of remote endpoints of the edges attached to v.
        // Not modified by the algorithm. 
        IntArrayList[] neighbend = new IntArrayList[nvertex];
        for (int i = 0; i < neighbend.length; i++) {
            neighbend[i] = new IntArrayList();
        }
        for (int k = 0; k < nedge; k++) {
            neighbend[edges[k][0]].add(2 * k + 1);
            neighbend[edges[k][1]].add(2 * k);
        }
        this.neighbend = neighbend;

        // If v is a vertex,
        // mate[v] is the remote endpoint of its matched edge, or -1 if it is single
        // (i.e. endpoint[mate[v]] is v's partner vertex).
        // Initially all vertices are single; updated during augmentation.
        int[] mate = new int[nvertex];
        Arrays.fill(mate, -1);
        this.mate = mate;

        // If b is a top-level blossom,
        // label[b] is 0 if b is unlabeled (free);
        //             1 if b is an S-vertex/blossom;
        //             2 if b is a T-vertex/blossom.
        // The label of a vertex is found by looking at the label of its
        // top-level containing blossom.
        // If v is a vertex inside a T-blossom,
        // label[v] is 2 iff v is reachable from an S-vertex outside the blossom.
        // Labels are assigned during a stage and reset after each augmentation.
        int[] label = new int[2 * nvertex];
        this.label = label;

        // If b is a labeled top-level blossom,
        // labelend[b] is the remote endpoint of the edge through which b obtained
        // its label, or -1 if b's base vertex is single.
        // If v is a vertex inside a T-blossom and label[v] == 2,
        // labelend[v] is the remote endpoint of the edge through which v is
        // reachable from outside the blossom.
        int[] labelend = new int[2 * nvertex];
        Arrays.fill(labelend, -1);
        this.labelend = labelend;

        // If v is a vertex,
        // inblossom[v] is the top-level blossom to which v belongs.
        // If v is a top-level vertex, v is itself a blossom (a trivial blossom)
        // and inblossom[v] == v.
        // Initially all vertices are top-level trivial blossoms.
        int[] inblossom = new int[nvertex];
        for (int i = 0; i < nvertex; i++) {
            inblossom[i] = i;
        }
        this.inblossom = inblossom;

        // If b is a sub-blossom,
        // blossomparent[b] is its immediate parent (sub-)blossom.
        // If b is a top-level blossom, blossomparent[b] is -1.
        int[] blossomparent = new int[2 * nvertex];
        Arrays.fill(blossomparent, -1);
        this.blossomparent = blossomparent;

        // If b is a non-trivial (sub-)blossom,
        // blossomchilds[b] is an ordered list of its sub-blossoms, starting with
        // the base and going round the blossom.
        IntArrayList[] blossomchilds = new IntArrayList[2 * nvertex];
        this.blossomchilds = blossomchilds;

        // If b is a (sub-)blossom,
        // blossombase[b] is its base VERTEX (i.e. recursive sub-blossom).
        int[] blossombase = new int[nvertex * 2];
        Arrays.fill(blossombase, -1);
        for (int i = 0; i < nvertex; i++) {
            blossombase[i] = i;
        }
        this.blossombase = blossombase;

        // If b is a non-trivial (sub-)blossom,
        // blossomendps[b] is a list of endpoints on its connecting edges,
        // such that blossomendps[b][i] is the local endpoint of blossomchilds[b][i]
        // on the edge that connects it to blossomchilds[b][wrap(i+1)].
        IntArrayList[] blossomendps = new IntArrayList[2 * nvertex];
        this.blossomendps = blossomendps;

        // If v is a free vertex (or an unreached vertex inside a T-blossom),
        // bestedge[v] is the edge to an S-vertex with least slack,
        // or -1 if there is no such edge.
        // If b is a (possibly trivial) top-level S-blossom,
        // bestedge[b] is the least-slack edge to a different S-blossom,
        // or -1 if there is no such edge.
        // This is used for efficient computation of delta2 and delta3.
        int[] bestedge = new int[2 * nvertex];
        Arrays.fill(bestedge, -1);
        this.bestedge = bestedge;

        // If b is a non-trivial top-level S-blossom,
        // blossombestedges[b] is a list of least-slack edges to neighbouring
        // S-blossoms, or None if no such list has been computed yet.
        // This is used for efficient computation of delta3.
        int[][] blossombestedges = new int[2 * nvertex][];
        this.blossombestedges = blossombestedges;

        // List of currently unused blossom numbers.
        IntArrayList unusedblossoms = new IntArrayList();
        for (int i = 0; i < nvertex; i++) {
            unusedblossoms.add(i + nvertex);
        }
        this.unusedblossoms = unusedblossoms;

        // If v is a vertex,
        // dualvar[v] = 2 * u(v) where u(v) is the v's variable in the dual
        // optimization problem (multiplication by two ensures integer values
        // throughout the algorithm if all edge weights are integers).
        // If b is a non-trivial blossom,
        // dualvar[b] = z(b) where z(b) is b's variable in the dual optimization
        // problem.
        int[] dualvar = new int[nvertex * 2];
        for (int i = 0; i < nvertex; i++) {
            dualvar[i] = maxweight;
        }
        this.dualvar = dualvar;

        // If allowedge[k] is true, edge k has zero slack in the optimization
        // problem; if allowedge[k] is false, the edge's slack may or may not
        // be zero.
        boolean[] allowedge = new boolean[nedge];
        this.allowedge = allowedge;

        // Queue of newly discovered S-vertices.
        IntArrayList queue = new IntArrayList();
        this.queue = queue;
    }

    // Return 2 * slack of edge k (does not work inside blossoms).
    private long slack(int k) {
        int[] edge = edges[k];
        return dualvar[edge[0]] + dualvar[edge[1]] - 2L * edge[2];
    }

    // Generate the leaf vertices of a blossom.
    private IntArrayList blossomLeaves(int b) {
        IntArrayList result = new IntArrayList();
        if (b < nvertex) {
            result.add(b);
        } else {
            for (int t : blossomchilds[b]) {
                if (t < nvertex) {
                    result.add(t);
                } else {
                    for (int v : blossomLeaves(t)) {
                        result.add(v);
                    }
                }
            }
        }
        return result;
    }

    // Assign label t to the top-level blossom containing vertex w
    // and record the fact that w was reached through the edge with
    // remote endpoint p.
    private void assignLabel(int w, int t, int p) {
        if (tracker != null) {
            tracker.start("func:assignLabel");
        }
        if (DEBUG) {
            System.out.println("assignLabel(" + w + "," + t + "," + p + ")");
        }
        int b = inblossom[w];
        assert label[w] == 0 && label[b] == 0;
        label[w] = label[b] = t;
        labelend[w] = labelend[b] = p;
        bestedge[w] = bestedge[b] = -1;
        if (t == 1) {
            // b became an S-vertex/blossom; add it(s vertices) to the queue.
            queue.addAll(blossomLeaves(b));
            if (DEBUG) {
                System.out.println("PUSH " + blossomLeaves(b));
            }
        } else if (t == 2) {
            // b became a T-vertex/blossom; assign label S to its mate.
            // (If b is a non-trivial blossom, its base is the only vertex
            // with an external mate.)
            int base = blossombase[b];
            assert mate[base] >= 0;
            assignLabel(endpoint[mate[base]], 1, mate[base] ^ 1);
        }
        if (tracker != null) {
            tracker.end();
        }
    }

    // Trace back from vertices v and w to discover either a new blossom
    // or an augmenting path. Return the base vertex of the new blossom or -1.
    private int scanBlossom(int v, int w) {
        if (DEBUG) {
            System.out.println("scanBlossom(" + v + "," + w + ")");
        }
        if (tracker != null) {
            tracker.start("func:scanBlossom");
        }
        // Trace back from v and w, placing breadcrumbs as we go.
        IntArrayList path = new IntArrayList();
        int base = -1;
        while (v != -1 || w != -1) {
            // Look for a breadcrumb in v's blossom or put a new breadcrumb.
            int b = inblossom[v];
            if ((label[b] & 4) != 0) {
                base = blossombase[b];
                break;
            }
            assert label[b] == 1;
            path.add(b);
            label[b] = 5;
            // Trace one step back.
            assert labelend[b] == mate[blossombase[b]];
            if (labelend[b] == -1) {
                // The base of blossom b is single; stop tracing this path.
                v = -1;
            } else {
                v = endpoint[labelend[b]];
                b = inblossom[v];
                assert label[b] == 2;
                // b is a T-blossom; trace one more step back.
                assert labelend[b] >= 0;
                v = endpoint[labelend[b]];
            }
            // Swap v and w so that we alternate between both paths.
            if (w != -1) {
                int t = v;
                v = w;
                w = t;
            }
        }
        // Remove breadcrumbs.
        for (int b : path) {
            label[b] = 1;
        }
        if (tracker != null) {
            tracker.end();
        }
        // Return base vertex, if we found one.
        return base;
    }

    // Construct a new blossom with given base, containing edge k which
    // connects a pair of S vertices. Label the new blossom as S; set its dual
    // variable to zero; relabel its T-vertices to S and add them to the queue.
    private void addBlossom(int base, int k) {
        if (tracker != null) {
            tracker.start("func:addBlossom");
        }
        int v = edges[k][0];
        int w = edges[k][1];
        int bb = inblossom[base];
        int bv = inblossom[v];
        int bw = inblossom[w];
        // Create blossom.
        int b = unusedblossoms.popInt();
        if (DEBUG) {
            System.out.println("addBlossom(" + base + "," + k + ") (v=" + v + " w=" + w + ") -> " + b);
        }
        blossombase[b] = base;
        blossomparent[b] = -1;
        blossomparent[bb] = b;
        // Make list of sub-blossoms and their interconnecting edge endpoints.
        IntArrayList path = new IntArrayList();
        blossomchilds[b] = path;
        IntArrayList endps = new IntArrayList();
        blossomendps[b] = endps;
        // Trace back from v to base.
        while (bv != bb) {
            // Add bv to the new blossom.
            blossomparent[bv] = b;
            path.add(bv);
            endps.add(labelend[bv]);
            assert (label[bv] == 2 || (label[bv] == 1 && labelend[bv] == mate[blossombase[bv]]));
            // Trace one step back.
            assert labelend[bv] >= 0;
            v = endpoint[labelend[bv]];
            bv = inblossom[v];
        }
        // Reverse lists, add endpoint that connects the pair of S vertices.
        path.add(bb);
        Collections.reverse(path);
        Collections.reverse(endps);
        endps.add(2 * k);
        // Trace back from w to base.
        while (bw != bb) {
            // Add bw to the new blossom.
            blossomparent[bw] = b;
            path.add(bw);
            endps.add(labelend[bw] ^ 1);
            assert (label[bw] == 2 || (label[bw] == 1 && labelend[bw] == mate[blossombase[bw]]));
            // Trace one step back.
            assert labelend[bw] >= 0;
            w = endpoint[labelend[bw]];
            bw = inblossom[w];
        }
        // Set label to S.
        assert label[bb] == 1;
        label[b] = 1;
        labelend[b] = labelend[bb];
        // Set dual variable to zero.
        dualvar[b] = 0;
        // Relabel vertices.
        for (int vi : blossomLeaves(b)) {
            if (label[inblossom[vi]] == 2) {
                // This T-vertex now turns into an S-vertex because it becomes
                // part of an S-blossom; add it to the queue.
                queue.add(vi);
            }
            inblossom[vi] = b;
        }
        // Compute blossombestedges[b].
        int[] bestedgeto = new int[2 * nvertex];
        Arrays.fill(bestedgeto, -1);
        for (int bvi : path) {
            int[][] nblists;
            if (blossombestedges[bvi] == null) {
                // This subblossom does not have a list of least-slack edges;
                // get the information from the vertices.
                IntArrayList blossomLeaves = blossomLeaves(bvi);
                nblists = new int[blossomLeaves.size()][];
                for (int i = 0; i < blossomLeaves.size(); i++) {
                    IntArrayList intArrayList = neighbend[blossomLeaves.getInt(i)];
                    nblists[i] = new int[intArrayList.size()];
                    for (int j = 0; j < intArrayList.size(); j++) {
                        nblists[i][j] = intArrayList.getInt(j) / 2;
                    }
                }
            } else {
                // Walk this subblossom's least-slack edges.
                nblists = new int[][]{blossombestedges[bvi]};
            }
            for (int[] nblist : nblists) {
                for (int ki : nblist) {
                    int i = edges[ki][0];
                    int j = edges[ki][1];
                    if (inblossom[j] == b) {
                        j = i;
                    }
                    int bj = inblossom[j];
                    if (bj != b && label[bj] == 1 && (bestedgeto[bj] == -1 || slack(ki) < slack(bestedgeto[bj]))) {
                        bestedgeto[bj] = ki;
                    }
                }
            }
            // Forget about least-slack edges of the subblossom.
            blossombestedges[bvi] = null;
            bestedge[bvi] = -1;
        }
        IntArrayList buffer = new IntArrayList();
        for (int ki : bestedgeto) {
            if (ki != -1) {
                buffer.add(ki);
            }
        }
        blossombestedges[b] = buffer.toIntArray();
        // Select bestedge[b].
        bestedge[b] = -1;
        for (int ki : blossombestedges[b]) {
            if (bestedge[b] == -1 || slack(ki) < slack(bestedge[b])) {
                bestedge[b] = ki;
            }
        }
        if (DEBUG) {
            System.out.println("blossomchilds[" + b + "]=" + blossomchilds[b]);
        }
        if (tracker != null) {
            tracker.end();
        }
    }

    // Expand the given top-level blossom.
    private void expandBlossom(int b, boolean endstage) {
        if (DEBUG) {
            System.out.println("expandBlossom(" + b + "," + endstage + ") " + blossomchilds[b]);
        }
        if (tracker != null) {
            tracker.start("func:expandBlossom");
        }
        // Convert sub-blossoms into top-level blossoms.
        for (int s : blossomchilds[b]) {
            blossomparent[s] = -1;
            if (s < nvertex) {
                inblossom[s] = s;
            } else if (endstage && dualvar[s] == 0) {
                // Recursively expand this sub-blossom.
                expandBlossom(s, endstage);
            } else {
                for (int v : blossomLeaves(s)) {
                    inblossom[v] = s;
                }
            }
        }
        // If we expand a T-blossom during a stage, its sub-blossoms must be
        // relabeled.
        if (!endstage && label[b] == 2) {
            // Start at the sub-blossom through which the expanding
            // blossom obtained its label, and relabel sub-blossoms untili
            // we reach the base.
            // Figure out through which sub-blossom the expanding blossom
            // obtained its label initially.
            assert labelend[b] >= 0;
            int entrychild = inblossom[endpoint[labelend[b] ^ 1]];
            // Decide in which direction we will go round the blossom.
            int j = blossomchilds[b].indexOf(entrychild);
            int endptrick;
            int jstep;
            if ((j & 1) != 0) {
                // Start index is odd; go forward and wrap.
                j -= blossomchilds[b].size();
                jstep = 1;
                endptrick = 0;
            } else {
                // Start index is even; go backward.
                jstep = -1;
                endptrick = 1;
            }
            // Move along the blossom until we get to the base.
            int p = labelend[b];
            while (j != 0) {
                // Relabel the T-sub-blossom.
                label[endpoint[p ^ 1]] = 0;
                label[endpoint[blossomendps[b].getInt((blossomendps[b].size() + (j - endptrick)) % blossomendps[b].size()) ^ endptrick ^ 1]] = 0;
                assignLabel(endpoint[p ^ 1], 2, p);
                // Step to the next S-sub-blossom and note its forward endpoint.
                allowedge[blossomendps[b].getInt((blossomendps[b].size() + (j - endptrick)) % blossomendps[b].size()) / 2] = true;
                j += jstep;
                p = blossomendps[b].getInt((blossomendps[b].size() + (j - endptrick)) % blossomendps[b].size()) ^ endptrick;
                // Step to the next T-sub-blossom.
                allowedge[p / 2] = true;
                j += jstep;
            }
            // Relabel the base T-sub-blossom WITHOUT stepping through to
            // its mate (so don't call assignLabel).
            int bv = blossomchilds[b].getInt(j);
            label[endpoint[p ^ 1]] = label[bv] = 2;
            labelend[endpoint[p ^ 1]] = labelend[bv] = p;
            bestedge[bv] = -1;
            // Continue along the blossom until we get back to entrychild.
            j += jstep;
            while (blossomchilds[b].getInt((blossomchilds[b].size() + j) % blossomchilds[b].size()) != entrychild) {
                // Examine the vertices of the sub-blossom to see whether
                // it is reachable from a neighbouring S-vertex outside the
                // expanding blossom.
                bv = blossomchilds[b].getInt((blossomchilds[b].size() + j) % blossomchilds[b].size());
                if (label[bv] == 1) {
                    // This sub-blossom just got label S through one of its
                    // neighbours; leave it.
                    j += jstep;
                    continue;
                }
                int v = -1;
                for (int vt : blossomLeaves(bv)) {
                    v = vt;
                    if (label[vt] != 0) {
                        break;
                    }
                }
                assert v != -1;
                // If the sub-blossom contains a reachable vertex, assign
                // label T to the sub-blossom.
                if (label[v] != 0) {
                    assert label[v] == 2;
                    assert inblossom[v] == bv;
                    label[v] = 0;
                    label[endpoint[mate[blossombase[bv]]]] = 0;
                    assignLabel(v, 2, labelend[v]);
                }
                j += jstep;
            }
        }
        // Recycle the blossom number.
        label[b] = labelend[b] = -1;
        blossomchilds[b] = blossomendps[b] = null;
        blossombase[b] = -1;
        blossombestedges[b] = null;
        bestedge[b] = -1;
        unusedblossoms.add(b);
        if (tracker != null) {
            tracker.end();
        }
    }

    // Swap matched/unmatched edges over an alternating path through blossom b
    // between vertex v and the base vertex. Keep blossom bookkeeping consistent.
    private void augmentBlossom(int b, int v) {
        if (DEBUG) {
            System.out.println("augmentBlossom(" + b + "," + v + ")");
        }
        if (tracker != null) {
            tracker.start("func:augmentBlossom");
        }
        // Bubble up through the blossom tree from vertex v to an immediate
        // sub-blossom of b.
        int t = v;
        while (blossomparent[t] != b) {
            t = blossomparent[t];
        }
        // Recursively deal with the first sub-blossom.
        if (t >= nvertex) {
            augmentBlossom(t, v);
        }
        // Decide in which direction we will go round the blossom.
        int i = blossomchilds[b].indexOf(t);
        int j = i;
        int jstep;
        int endptrick;
        if ((i & 1) != 0) {
            // Start index is odd; go forward and wrap.
            j -= blossomchilds[b].size();
            jstep = 1;
            endptrick = 0;
        } else {
            // Start index is even; go backward.
            jstep = -1;
            endptrick = 1;
        }
        // Move along the blossom until we get to the base.
        while (j != 0) {
            // Step to the next sub-blossom and augment it recursively.
            j += jstep;
            t = blossomchilds[b].getInt((blossomchilds[b].size() + j) % blossomchilds[b].size());
            int p = blossomendps[b].getInt((blossomendps[b].size() + (j - endptrick)) % blossomendps[b].size()) ^ endptrick;
            if (t >= nvertex) {
                augmentBlossom(t, endpoint[p]);
            }
            // Step to the next sub-blossom and augment it recursively.
            j += jstep;
            t = blossomchilds[b].getInt((blossomchilds[b].size() + j) % blossomchilds[b].size());
            if (t >= nvertex) {
                augmentBlossom(t, endpoint[p ^ 1]);
            }
            // Match the edge connecting those sub-blossoms.
            mate[endpoint[p]] = p ^ 1;
            mate[endpoint[p ^ 1]] = p;
            if (DEBUG) {
                System.out.println("PAIR " + endpoint[p] + " " + endpoint[p ^ 1] + " (k=" + p / 2 + ")");
            }
        }
        // Rotate the list of sub-blossoms to put the new base at the front.
        for (int c = 0; c < i; c++) {
            blossomchilds[b].add(blossomchilds[b].getInt(c));
            blossomendps[b].add(blossomendps[b].getInt(c));
        }
        blossomchilds[b].removeElements(0, i);
        blossomendps[b].removeElements(0, i);
        blossombase[b] = blossombase[blossomchilds[b].getInt(0)];
        assert blossombase[b] == v;
        if (tracker != null) {
            tracker.end();
        }
    }

    // Swap matched/unmatched edges over an alternating path between two
    // single vertices. The augmenting path runs through edge k, which
    // connects a pair of S vertices.
    private void augmentMatching(int k) {
        if (tracker != null) {
            tracker.start("func:augmentMatching");
        }
        int v = edges[k][0];
        int w = edges[k][1];
        if (DEBUG) {
            System.out.println("augmentMatching(" + k + ") (v=" + v + " w=" + w + ")");
            System.out.println("PAIR " + v + " " + w + " (k=" + k + ")");
        }
        for (int[] sp : new int[][]{new int[]{v, 2 * k + 1}, new int[]{w, 2 * k}}) {
            int s = sp[0];
            int p = sp[1];
            // Match vertex s to remote endpoint p. Then trace back from s
            // until we find a single vertex, swapping matched and unmatched
            // edges as we go.
            while (true) {
                int bs = inblossom[s];
                assert label[bs] == 1;
                assert labelend[bs] == mate[blossombase[bs]];
                // Augment through the S-blossom from s to base.
                if (bs >= nvertex) {
                    augmentBlossom(bs, s);
                }
                // Update mate[s]
                mate[s] = p;
                // Trace one step back.
                if (labelend[bs] == -1) {
                    // Reached single vertex; stop.
                    break;
                }
                int t = endpoint[labelend[bs]];
                int bt = inblossom[t];
                assert label[bt] == 2;
                // Trace one step back.
                assert labelend[bt] >= 0;
                s = endpoint[labelend[bt]];
                int j = endpoint[labelend[bt] ^ 1];
                // Augment through the T-blossom from j to base.
                assert blossombase[bt] == t;
                if (bt >= nvertex) {
                    augmentBlossom(bt, j);
                }
                // Update mate[j]
                mate[j] = labelend[bt];
                // Keep the opposite endpoint;
                // it will be assigned to mate[s] in the next step.
                p = labelend[bt] ^ 1;
                if (DEBUG) {
                    System.out.println("PAIR " + s + " " + t + "(k=" + (p / 2) + ")");
                }
            }
        }
        if (tracker != null) {
            tracker.end();
        }
    }

    // Verify that the optimum solution has been reached.
    private void verifyOptimum() {
        int vdualoffset;
        int dualvarminleft = Integer.MAX_VALUE;
        for (int i = 0; i < nvertex; i++) {
            dualvarminleft = Math.min(dualvarminleft, dualvar[i]);
        }
        int dualvarminright = Integer.MAX_VALUE;
        for (int i = nvertex; i < dualvar.length; i++) {
            dualvarminright = Math.min(dualvarminright, dualvar[i]);
        }
        if (maxcardinality) {
            // Vertices may have negative dual;
            // find a constant non-negative number to add to all vertex duals.
            vdualoffset = Math.max(0, -dualvarminleft);
        } else {
            vdualoffset = 0;
        }
        // 0. all dual variables are non-negative
        assert dualvarminleft + vdualoffset >= 0;
        assert dualvarminright >= 0;
        // 0. all edges have non-negative slack and
        // 1. all matched edges have zero slack;
        for (int k = 0; k < nedge; k++) {
            int i = edges[k][0];
            int j = edges[k][1];
            int wt = edges[k][2];
            long s = dualvar[i] + dualvar[j] - 2L * wt;
            IntArrayList iblossoms = new IntArrayList(new int[]{i});
            IntArrayList jblossoms = new IntArrayList(new int[]{j});
            while (blossomparent[iblossoms.getInt(iblossoms.size() - 1)] != -1) {
                iblossoms.add(blossomparent[iblossoms.getInt(iblossoms.size() - 1)]);
            }
            while (blossomparent[jblossoms.getInt(jblossoms.size() - 1)] != -1) {
                jblossoms.add(blossomparent[jblossoms.getInt(jblossoms.size() - 1)]);
            }
            Collections.reverse(iblossoms);
            Collections.reverse(jblossoms);
            assert iblossoms.size() == jblossoms.size();
            for (int c = 0; c < iblossoms.size(); c++) {
                int bi = iblossoms.getInt(c);
                int bj = jblossoms.getInt(c);
                if (bi != bj) {
                    break;
                }
                s += 2 * dualvar[bi];
            }
            assert s >= 0;
            if (mate[i] / 2 == k || mate[j] / 2 == k) {
                assert mate[i] / 2 == k && mate[j] / 2 == k;
                assert s == 0;
            }
        }
        // 2. all single vertices have zero dual value;
        for (int v = 0; v < nvertex; v++) {
            assert mate[v] >= 0 || dualvar[v] + vdualoffset == 0;
        }
        // 3. all blossoms with positive dual value are full.
        for (int b = nvertex; b < 2 * nvertex; b++) {
            if (blossombase[b] >= 0 && dualvar[b] > 0) {
                assert blossomendps[b].size() % 2 == 1;
                for (int i = 1; i < blossomendps[b].size(); i += 2) {
                    int p = blossomendps[b].getInt(i);
                    assert mate[endpoint[p]] == (p ^ 1);
                    assert mate[endpoint[p ^ 1]] == p;
                }
            }
        }
        // Ok.
    }

    // Check optimized delta2 against a trivial computation.
    private void checkDelta2() {
        for (int v = 0; v < nvertex; v++) {
            if (label[inblossom[v]] == 0) {
                long bd = Long.MAX_VALUE;
                int bk = -1;
                for (int p : neighbend[v]) {
                    int k = p / 2;
                    int w = endpoint[p];
                    if (label[inblossom[w]] == 1) {
                        long d = slack(k);
                        if (bk == -1 || d < bd) {
                            bk = k;
                            bd = d;
                        }
                    }
                }
                if (DEBUG) {
                    if ((bestedge[v] != -1 || bk != -1) && (bestedge[v] == -1 || bd != slack(bestedge[v]))) {
                        System.out.println(
                                "v=" + v + " bk=" + bk + " bd=" + bd + " bestedge=" + bestedge[v] + " slack=" + slack(bestedge[v])
                        );
                    }
                }
                assert (bk == -1 && bestedge[v] == -1) || (bestedge[v] != -1 && bd == slack(bestedge[v]));
            }
        }
    }

    // Check optimized delta3 against a trivial computation.
    private void checkDelta3() {
        int bk = -1;
        long bd = Long.MAX_VALUE;
        int tbk = -1;
        long tbd = Long.MAX_VALUE;
        for (int b = 0; b < 2 * nvertex; b++) {
            if (blossomparent[b] == -1 && label[b] == 1) {
                for (int v : blossomLeaves(b)) {
                    for (int p : neighbend[v]) {
                        int k = p / 2;
                        int w = endpoint[p];
                        if (inblossom[w] != b && label[inblossom[w]] == 1) {
                            long d = slack(k);
                            if (bk == -1 || d < bd) {
                                bk = k;
                                bd = d;
                            }
                        }
                    }
                }
                if (bestedge[b] != -1) {
                    int i = edges[bestedge[b]][0];
                    int j = edges[bestedge[b]][1];
                    assert inblossom[i] == b || inblossom[j] == b;
                    assert inblossom[i] != b || inblossom[j] != b;
                    assert label[inblossom[i]] == 1 && label[inblossom[j]] == 1;
                    if (tbk == -1 || slack(bestedge[b]) < tbd) {
                        tbk = bestedge[b];
                        tbd = slack(bestedge[b]);
                    }
                }
            }
        }
        if (DEBUG) {
            if (bd != tbd) {
                System.out.println("bk=" + bk + " tbk=" + tbk + " bd=" + bd + " tbd=" + tbd);
            }
        }
        assert bd == tbd;
    }

    /* Compute a maximum-weighted matching in the general undirected
    weighted graph given by "edges".  If "maxcardinality" is true,
    only maximum-cardinality matchings are considered as solutions.

    Edges is a sequence of tuples (i, j, wt) describing an undirected
    edge between vertex i and vertex j with weight wt.  There is at most
    one edge between any two vertices; no vertex has an edge to itself.
    Vertices are identified by consecutive, non-negative integers.

    Return a list "mate", such that mate[i] == j if vertex i is
    matched to vertex j, and mate[i] == -1 if vertex i is not matched.

    This function takes time O(n ** 3). */
    public int[] maxWeightMatching() {
        // Deal swiftly with empty graphs.
        if (this.edges.length == 0) {
            return new int[0];
        }

        // Main loop: continue until no further improvement is possible.
        for (int t = 0; t < nvertex; t++) {

            // Each iteration of this loop is a "stage".
            // A stage finds an augmenting path and uses that to improve
            // the matching.
            if (DEBUG) {
                System.out.println("STAGE " + t);
            }

            // Remove labels from top-level blossoms/vertices.
            Arrays.fill(label, 0);

            // Forget all about least-slack edges.
            Arrays.fill(bestedge, -1);
            for (int i = nvertex; i < blossombestedges.length; i++) {
                blossombestedges[i] = null;
            }

            // Loss of labeling means that we can not be sure that currently
            // allowable edges remain allowable througout this stage.
            Arrays.fill(allowedge, false);

            // Make queue empty.
            queue.clear();

            // Label single blossoms/vertices with S and put them in the queue.
            for (int v = 0; v < nvertex; v++) {
                if (mate[v] == -1 && label[inblossom[v]] == 0) {
                    assignLabel(v, 1, -1);
                }
            }

            // Loop until we succeed in augmenting the matching.
            int augmented = 0;
            while (true) {
                // Each iteration of this loop is a "substage".
                // A substage tries to find an augmenting path;
                // if found, the path is used to improve the matching and
                // the stage ends. If there is no augmenting path, the
                // primal-dual method is used to pump some slack out of
                // the dual variables.
                if (DEBUG) {
                    System.out.println("SUBSTAGE");
                }

                if (tracker != null) {
                    tracker.start("labelling");
                }
                // Continue labeling until all vertices which are reachable
                // through an alternating path have got a label.
                while (!queue.isEmpty() && augmented == 0) {

                    // Take an S vertex from the queue.
                    int v = queue.popInt();
                    if (DEBUG) {
                        System.out.println("POP v=" + v);
                    }
                    assert label[inblossom[v]] == 1;

                    // Scan its neighbours:
                    int w;
                    for (int p : neighbend[v]) {
                        int k = p / 2;
                        w = endpoint[p];
                        // w is a neighbour to v
                        if (inblossom[v] == inblossom[w]) {
                            // this edge is internal to a blossom; ignore it
                            continue;
                        }
                        long kslack = Long.MAX_VALUE;
                        if (!allowedge[k]) {
                            kslack = slack(k);
                            if (kslack <= 0) {
                                // edge k has zero slack => it is allowable
                                allowedge[k] = true;
                            }
                        }
                        if (allowedge[k]) {
                            if (label[inblossom[w]] == 0) {
                                // (C1) w is a free vertex;
                                // label w with T and label its mate with S (R12).
                                assignLabel(w, 2, p ^ 1);
                            } else if (label[inblossom[w]] == 1) {
                                // (C2) w is an S-vertex (not in the same blossom);
                                // follow back-links to discover either an
                                // augmenting path or a new blossom.
                                int base = scanBlossom(v, w);
                                if (base >= 0) {
                                    // Found a new blossom; add it to the blossom
                                    // bookkeeping and turn it into an S-blossom.
                                    addBlossom(base, k);
                                } else {
                                    // Found an augmenting path; augment the
                                    // matching and end this stage.
                                    augmentMatching(k);
                                    augmented = 1;
                                    break;
                                }
                            } else if (label[w] == 0) {
                                // w is inside a T-blossom, but w itself has not
                                // yet been reached from outside the blossom;
                                // mark it as reached (we need this to relabel
                                // during T-blossom expansion).
                                assert label[inblossom[w]] == 2;
                                label[w] = 2;
                                labelend[w] = p ^ 1;
                            }
                        } else if (label[inblossom[w]] == 1) {
                            // keep track of the least-slack non-allowable edge to
                            // a different S-blossom.
                            int b = inblossom[v];
                            if (bestedge[b] == -1 || kslack < slack(bestedge[b])) {
                                bestedge[b] = k;
                            }
                        } else if (label[w] == 0) {
                            // w is a free vertex (or an unreached vertex inside
                            // a T-blossom) but we can not reach it yet;
                            // keep track of the least-slack edge that reaches w.
                            if (bestedge[w] == -1 || kslack < slack(bestedge[w])) {
                                bestedge[w] = k;
                            }
                        }
                    }
                }
                if (tracker != null) {
                    tracker.end();
                }

                if (augmented != 0) {
                    break;
                }

                // There is no augmenting path under these constraints;
                // compute delta and reduce slack in the optimization problem.
                // (Note that our vertex dual variables, edge slacks and delta's
                // are pre-multiplied by two.)
                int deltatype = -1;
                long delta = Long.MAX_VALUE;
                int deltaedge = Integer.MAX_VALUE;
                int deltablossom = Integer.MAX_VALUE;

                // Verify data structures for delta2/delta3 computation.
                if (CHECK_DELTA) {
                    checkDelta2();
                    checkDelta3();
                }

                // Compute delta1: the minumum value of any vertex dual.
                if (!maxcardinality) {
                    deltatype = 1;
                    for (int i = 0; i < nvertex; i++) {
                        delta = Math.min(delta, dualvar[i]);
                    }
                }

                if (tracker != null) {
                    tracker.start("delta2");
                }
                // Compute delta2: the minimum slack on any edge between
                // an S-vertex and a free vertex.
                for (int v = 0; v < nvertex; v++) {
                    if (label[inblossom[v]] == 0 && bestedge[v] != -1) {
                        long d = slack(bestedge[v]);
                        if (deltatype == -1 || d < delta) {
                            delta = d;
                            deltatype = 2;
                            deltaedge = bestedge[v];
                        }
                    }
                }
                if (tracker != null) {
                    tracker.end();
                }

                if (tracker != null) {
                    tracker.start("delta3");
                }
                // Compute delta3: half the minimum slack on any edge between
                // a pair of S-blossoms.
                for (int b = 0; b < 2 * nvertex; b++) {
                    if (blossomparent[b] == -1 && label[b] == 1 && bestedge[b] != -1) {
                        long kslack = slack(bestedge[b]);
                        long d = kslack / 2;
                        if (deltatype == -1 || d < delta) {
                            delta = d;
                            deltatype = 3;
                            deltaedge = bestedge[b];
                        }
                    }
                }
                if (tracker != null) {
                    tracker.end();
                }

                if (tracker != null) {
                    tracker.start("delta4");
                }
                // Compute delta4: minimum z variable of any T-blossom.
                for (int b = nvertex; b < 2 * nvertex; b++) {
                    if (blossombase[b] >= 0 && blossomparent[b] == -1 && label[b] == 2 && (deltatype == -1 || dualvar[b] < delta)) {
                        delta = dualvar[b];
                        deltatype = 4;
                        deltablossom = b;
                    }
                }
                if (tracker != null) {
                    tracker.end();
                }

                if (deltatype == -1) {
                    // No further improvement possible; max-cardinality optimum
                    // reached. Do a final delta update to make the optimum
                    // verifyable.
                    assert maxcardinality;
                    deltatype = 1;
                    int mindualvar = Integer.MAX_VALUE;
                    for (int i = 0; i < nvertex; i++) {
                        mindualvar = Math.min(mindualvar, dualvar[i]);
                    }
                    delta = Math.max(0, mindualvar);
                }

                if (tracker != null) {
                    tracker.start("update_dual");
                }
                // Update dual variables according to delta.
                for (int v = 0; v < nvertex; v++) {
                    if (label[inblossom[v]] == 1) {
                        // S-vertex: 2*u = 2*u - 2*delta
                        dualvar[v] -= delta;
                    } else if (label[inblossom[v]] == 2) {
                        // T-vertex: 2*u = 2*u + 2*delta
                        dualvar[v] += delta;
                    }
                }
                for (int b = nvertex; b < 2 * nvertex; b++) {
                    if (blossombase[b] >= 0 && blossomparent[b] == -1) {
                        if (label[b] == 1) {
                            // top-level S-blossom: z = z + 2*delta
                            dualvar[b] += delta;
                        } else if (label[b] == 2) {
                            // top-level T-blossom: z = z - 2*delta
                            dualvar[b] -= delta;
                        }
                    }
                }
                if (tracker != null) {
                    tracker.end();
                }

                // Take action at the point where minimum delta occurred.
                if (DEBUG) {
                    System.out.println("delta" + deltatype + "=" + delta);
                }
                int i;
                int j;
                if (deltatype == 1) {
                    // No further improvement possible; optimum reached.
                    break;
                } else if (deltatype == 2) {
                    // Use the least-slack edge to continue the search.
                    allowedge[deltaedge] = true;
                    i = edges[deltaedge][0];
                    j = edges[deltaedge][1];
                    if (label[inblossom[i]] == 0) {
                        int ti = i;
                        i = j;
                        j = ti;
                    }
                    assert label[inblossom[i]] == 1;
                    queue.add(i);
                } else if (deltatype == 3) {
                    // Use the least-slack edge to continue the search.
                    allowedge[deltaedge] = true;
                    i = edges[deltaedge][0];
                    j = edges[deltaedge][1];
                    assert label[inblossom[i]] == 1;
                    queue.add(i);
                } else if (deltatype == 4) {
                    // Expand the least-z blossom.
                    expandBlossom(deltablossom, false);
                }

                // End of a this substage.
            }

            // Stop when no more augmenting path can be found.
            if (augmented == 0) {
                break;
            }

            // End of a stage; expand all S-blossoms which have dualvar = 0.
            for (int b = nvertex; b < 2 * nvertex; b++) {
                if (blossomparent[b] == -1 && blossombase[b] >= 0 && label[b] == 1 && dualvar[b] == 0) {
                    expandBlossom(b, true);
                }
            }
        }
        // Verify that we reached the optimum solution.
        if (CHECK_OPTIMUM) {
            verifyOptimum();
        }

        // Transform mate[] such that mate[v] is the vertex to which v is paired.
        for (int v = 0; v < nvertex; v++) {
            if (mate[v] >= 0) {
                mate[v] = endpoint[mate[v]];
            }
        }
        for (int v = 0; v < nvertex; v++) {
            assert mate[v] == -1 || mate[mate[v]] == v;
        }

        if (tracker != null) {
            System.out.println(tracker);
        }

        return mate;
    }

}
