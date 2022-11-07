
public class FakeClusteringPipeline {
    
    public static void main(String[] args) {
        String outString = "java -jar clustering.jar";
        for (String arg: args) {
            outString += " " + arg;
        }
        System.out.println(outString);
    }

}
