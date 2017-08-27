import it.unimi.dsi.fastutil.longs.LongArrayList;
import it.unimi.dsi.fastutil.objects.Object2ObjectMap;
import it.unimi.dsi.fastutil.objects.Object2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;

import java.util.Comparator;

public class Tracker {

    private ObjectArrayList<String> funcStack = new ObjectArrayList<>();
    private LongArrayList timeStack = new LongArrayList();
    private Object2ObjectOpenHashMap<String, int[]> stats = new Object2ObjectOpenHashMap<>();

    public final void start(String func) {
        if (!stats.containsKey(func)) {
            stats.put(func, new int[2]);
        }
        funcStack.add(func);
        stats.get(func)[0]++;
        timeStack.add(System.nanoTime());
    }

    public final void end() {
        stats.get(funcStack.pop())[1] += System.nanoTime() - timeStack.popLong();
    }

    @Override
    public final String toString() {
        assert timeStack.isEmpty();
        assert funcStack.isEmpty();
        ObjectArrayList<Object2ObjectMap.Entry<String, int[]>> entryList = new ObjectArrayList<>();
        entryList.addAll(stats.object2ObjectEntrySet());
        entryList.sort(Comparator.comparingInt(o -> o.getValue()[0]));

        StringBuilder result = new StringBuilder();
        for (Object2ObjectMap.Entry<String, int[]> entry : entryList) {
            result.append(entry.getKey()).append(" : ").append(entry.getValue()[0]).append(" ").append(entry.getValue()[1] / 1000000.0).append("\n");
        }
        return result.toString();
    }

}
