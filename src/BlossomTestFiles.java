import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import org.junit.Assert;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestName;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.stream.Collectors;


@RunWith(Parameterized.class)
public class BlossomTestFiles {

    private static final String testFolder = "tests";

    @Rule
    public TestName name = new TestName();

    private static ObjectArrayList<Object[]> getTestFileNames() {
        return Arrays.stream(
                new File(testFolder).listFiles()
        ).filter(
                File::isFile
        ).map(
                File::getName
        ).filter(
                x -> x.startsWith("in_")
        ).distinct().sorted().map(
                fileName -> new String[]{fileName}
        ).collect(
                Collectors.toCollection(ObjectArrayList::new)
        );
    }

    private String getTestName() {
        return name.getMethodName().replaceAll("(.*?\\[)?([^\\[]*?)\\]?", "$2");
    }

    @Parameterized.Parameters(name = "{0}")
    public static Collection<Object[]> data() {
        return getTestFileNames();
    }

    public BlossomTestFiles(String ignored) { }

    @Test
    public void test() throws IOException {
        String testName = getTestName();
        ObjectArrayList<int[]> edgesList = new ObjectArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(new File(testFolder, testName)))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] entries = line.split(" ");
                edgesList.add(new int[] {Integer.valueOf(entries[0]), Integer.valueOf(entries[1]), Integer.valueOf(entries[2])});
            }
        }
        int[][] edges = new int[edgesList.size()][];
        edgesList.toArray(edges);
        int[] matches = new Blossom(edges, true).maxWeightMatching();
        int[] expectedMatches = new int[matches.length];
        try (BufferedReader br = new BufferedReader(new FileReader(new File(testFolder, testName.replace("in_", "out_"))))) {
            String line;
            int i = 0;
            while ((line = br.readLine()) != null) {
                expectedMatches[i++] = Integer.valueOf(line);
            }
        }
        Assert.assertArrayEquals(matches, expectedMatches);
    }

}