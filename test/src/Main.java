
import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.StringTokenizer;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Enumeration;

public class Main {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		//RandomAccessFile goldstandardFile = OpenFile("/n/staffstore/burcu/POSCORPUS_tags");
		RandomAccessFile goldstandardFile = OpenFile("pennTagged500_17.txt");
		RandomAccessFile resultsFile = OpenFile("taggedp.txt");
		String line;
		
		int corpusSize=0;
		ArrayList<String> words = new ArrayList<String>();
		
		while((line = ReadLine(resultsFile)) != null){
			int index = line.lastIndexOf(" ");
//			tag = "t" + (Integer.parseInt(wordd.substring(index+1))+1);
			line = line.substring(0,index);
			//line=line.split(" ",2)[1];
//			System.out.println(line);
			StringTokenizer tokenizer = new StringTokenizer(line.toLowerCase(), " ");
			
			while(tokenizer.hasMoreTokens()){
				String nextToken = tokenizer.nextToken();
//				System.out.println(nextToken);
				if(!nextToken.equals("SB_BURCU/-1")){
					words.add(nextToken);
					corpusSize++;
				}
			}
		}
		System.out.println("Corpus size:" + corpusSize);
		
		ArrayList<String> tags = new ArrayList<String>();
		int counter = 0;
		while((line = ReadLine(goldstandardFile)) != null && counter<=corpusSize){
			StringTokenizer tokenizer = new StringTokenizer(line.toLowerCase(), " ");
			
			while(tokenizer.hasMoreTokens() && counter<=corpusSize){
				String nextToken = tokenizer.nextToken();
				tags.add(nextToken);
				counter++;
			}
		}

		Hashtable<String, String> idmatches = MatchTags(words, tags);
		
		double accuracy = Evaluate(idmatches, words, tags);
		
		System.out.println("Accuracy: " + accuracy);
		
		CloseFile(goldstandardFile);
		CloseFile(resultsFile);
	}
	
	public static double Evaluate(Hashtable<String, String> idmatches, ArrayList<String> words, ArrayList<String> tags){
		
		int correct = 0;
		
		for(int i=0; i<words.size(); i++){
			
			String nextResultWord	= words.get(i);
			String nextGoldWord		= tags.get(i);

			StringTokenizer tokenizer = new StringTokenizer(nextResultWord, "/");
			String resultWord		= tokenizer.nextToken();
			String resultPosId	= tokenizer.nextToken();
			if(tokenizer.hasMoreTokens())
				resultPosId = tokenizer.nextToken();
			
			tokenizer = new StringTokenizer(nextGoldWord, "/");
			String goldWord		= tokenizer.nextToken();
			String goldPosId	= tokenizer.nextToken();
	
			StringTokenizer tokenizer3 = new StringTokenizer(goldPosId, "|");
			while(tokenizer3.hasMoreTokens()){
				String next = tokenizer3.nextToken();
				if(idmatches.get(resultPosId).equals(next))
					correct++;
			}
		}
		return (double)100*correct/(double)(words.size());
	}
	
	public static Hashtable<String, String> MatchTags(ArrayList<String> words, ArrayList<String> tags){
		
		Hashtable<String, Hashtable<String, Integer>> idmatches = new Hashtable<String, Hashtable<String, Integer>>();  
		
		for(int i=0; i<words.size()-1; i++){
			
			String nextResultWord	= words.get(i);
			String nextGoldWord		= tags.get(i);
			//System.out.println(nextResultWord);
			StringTokenizer tokenizer = new StringTokenizer(nextResultWord, "/");
			String resultWord		= tokenizer.nextToken();
			String resultPosId	= tokenizer.nextToken();
			//System.out.println(nextResultWord);
			if(tokenizer.hasMoreTokens())
				resultPosId = tokenizer.nextToken();
			
			tokenizer = new StringTokenizer(nextGoldWord, "/");
			String goldWord		= tokenizer.nextToken();
			String goldPosId	= tokenizer.nextToken();
			
			ArrayList<String> goldposs = new ArrayList<String>();
			StringTokenizer tokenizer3 = new StringTokenizer(goldPosId, "|");
			while(tokenizer3.hasMoreTokens()){
				String next = tokenizer3.nextToken();
				goldposs.add(next);
			}
			
			if(idmatches.containsKey(resultPosId)){
				
				Hashtable<String, Integer> realids = idmatches.get(resultPosId);
				
				for(int g=0; g<goldposs.size(); g++){
					if(realids.containsKey(goldposs.get(g))){
						int freq = realids.get(goldposs.get(g));
						realids.put(goldposs.get(g), freq+1);
					}
					else
						realids.put(goldposs.get(g), 1);
				}
			}
			else{
				Hashtable<String, Integer> newrealids = new Hashtable<String, Integer>();
				
				for(int g=0; g<goldposs.size(); g++)
					newrealids.put(goldposs.get(g), 1);
				
				idmatches.put(resultPosId, newrealids);
			}
		}
		
		Hashtable<String, String> matches = new Hashtable<String, String>(); 
		
		Enumeration<String> k = idmatches.keys();
		while (k.hasMoreElements()) {
			String key = (String) k.nextElement();
			
			int max = 0;
			String maxId = null;
			Hashtable<String, Integer> realids = idmatches.get(key);
			
			Enumeration<String> rid = realids.keys();
			while (rid.hasMoreElements()) {
				String key2 = rid.nextElement();
				int freq = realids.get(key2);
				if(freq>max){
					maxId = key2;
					max = freq;
				}
			}
			matches.put(key, maxId);
		}
		
		return matches;
	}
	
	public static RandomAccessFile OpenFile(String fileName) {
		// TODO Auto-generated method stub
		try{
			File file = new File(fileName);
			RandomAccessFile raf = new RandomAccessFile(file, "rw");
			return raf;
			}
		
		catch (IOException e) {
			System.out.println("Unable to open "+ fileName +": " + e.getMessage());
			System.out.println("IOException:");
			e.printStackTrace();
	        }
		return null;
	}
	
	public static String ReadLine(RandomAccessFile raf) {
		// TODO Auto-generated method stub
		try{
			return raf.readLine();
		}
		catch (IOException e) {
			System.out.println("IOException:");
			e.printStackTrace();
		}
		return null;
	}
	
	public static void CloseFile(RandomAccessFile raf) {
		// TODO Auto-generated method stub
		try{
			raf.close();
			
		}
		catch (IOException e) {
			System.out.println("IOException:");
			e.printStackTrace();
		}
	}
	
	public static void WriteString(String str, RandomAccessFile raf) {
		// TODO Auto-generated method stub
		try{
			raf.writeBytes(str);
		}
		catch (IOException e) {
			System.out.println("IOException:");
			e.printStackTrace();
		}
	}
}
