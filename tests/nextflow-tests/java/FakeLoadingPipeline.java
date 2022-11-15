
public class FakeLoadingPipeline {
    
    public static void main(String[] args) {
        String outString = "java -jar loading.jar";
        for (String arg: args) {
            outString += " " + arg;
        }
        System.out.println(outString);
    }

}
