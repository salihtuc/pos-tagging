import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Properties;

public class MainRead {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Properties prop = new Properties();
		InputStream input = null;
		String paramsFile;
		String inputFile;
		int tagsize=0;
		if (args.length > 0) {
			paramsFile = args[0];
		} else {
			paramsFile = "params.properties";
		}
		// get params
		try {

			input = new FileInputStream(paramsFile);

			// load a properties file
			prop.load(input);
			inputFile=prop.getProperty("inputFile");
			tagsize=Integer.parseInt(prop.getProperty("tagSize"));
			
			readText(inputFile);
			
			System.out.println(sentences);

		} catch (IOException ex) {
			ex.printStackTrace();
		} finally {
			if (input != null) {
				try {
					input.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	static ArrayList<String> sentences=new ArrayList<String>();
	public static void readText(String inputFile) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(inputFile));
		
		 try {
	            StringBuilder sb;
	            String line = br.readLine();
	            while (line != null) {
	            	sentences.add(line.toLowerCase());
	                line = br.readLine();
	            }
	        } finally {
	            br.close();
	        }
	        System.out.println("Read in "+sentences.size()+" sentences");
	}

}
