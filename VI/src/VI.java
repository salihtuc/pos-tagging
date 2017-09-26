import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.StringTokenizer;


public class VI {

	static Hashtable<String, Hashtable<String, Integer>> resultClusters = new Hashtable<String, Hashtable<String, Integer>>();
	static Hashtable<String, Hashtable<String, Integer>> goldClusters = new Hashtable<String, Hashtable<String, Integer>>();
	static Hashtable<String, Integer> resultClustersSize = new Hashtable<String, Integer>();
	static Hashtable<String, Integer> goldClustersSize = new Hashtable<String, Integer>();
	static int corpusSize=0;
	static int counter = 0;
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		

		
		//RandomAccessFile goldstandardFile = OpenFile("/n/staffstore/burcu/POSCORPUS_tags");
		RandomAccessFile goldstandardFile = OpenFile("tagged400.txt");
		RandomAccessFile resultsFile = OpenFile("deneme_out400.txt");
		String line=null, word=null;
		//Hashtable<String, String> clusters = new Hashtable<String, String>();
		//Hashtable<String, String> postags = new Hashtable<String, String>();
		Hashtable<String, String> wordss = new Hashtable<String, String>();
		
		while((line = ReadLine(resultsFile)) != null){
			StringTokenizer tokenizer = new StringTokenizer(line, " ");
			
			while(tokenizer.hasMoreTokens()){
				String nextToken = tokenizer.nextToken();
				
				if(!nextToken.equals("SB_BURCU/-1")){
					StringTokenizer tokenizer2 = new StringTokenizer(nextToken, "/");
					word = tokenizer2.nextToken();
//					System.out.println(nextToken);
					String nextCluster = tokenizer2.nextToken();
//					System.out.println(nextToken);
					if(tokenizer2.hasMoreTokens()){
						word += nextCluster;
						nextCluster = tokenizer2.nextToken();
					}
					wordss.put(word, word);
					AddToHashHash(nextCluster, word, resultClusters);
					AddToHashtable_str(resultClustersSize, nextCluster);
					corpusSize++;
				}
			}
		}
		
		//while((line = ReadLine(goldstandardFile)) != null && counter<=corpusSize){
		while((line = ReadLine(goldstandardFile)) != null && counter <= corpusSize){
			StringTokenizer tokenizer = new StringTokenizer(line, " ");
			
			while(tokenizer.hasMoreTokens() && counter <= corpusSize){
				String nextToken = tokenizer.nextToken();
				StringTokenizer tokenizer2 = new StringTokenizer(nextToken, "/");
				word = tokenizer2.nextToken();
				String nextpostag = tokenizer2.nextToken();
				
				StringTokenizer tokenizer3 = new StringTokenizer(nextpostag, "|");
				while(tokenizer3.hasMoreTokens() && counter <= corpusSize){
					String next = tokenizer3.nextToken();
					AddToHashHash(next, word, goldClusters);
					AddToHashtable_str(goldClustersSize, next);
					counter++;
					if(tokenizer3.hasMoreTokens())
						tokenizer3.nextToken();
				}
				
			//	postags.put(nextpostag, nextpostag);   
			}
		}
		
		double resultsEntropy = CalEntropy(resultClusters);
		double goldEntropy = CalEntropyGold(goldClusters);
		
		double calMutual = CalMutualInformation();
		
		System.out.println("VI "+(resultsEntropy+goldEntropy+(-2)*calMutual));
		
		System.out.println("NMI "+calMutual/Math.sqrt(resultsEntropy*goldEntropy));
		
		double completeness = 1-(resultsEntropy-calMutual)/resultsEntropy;
		double homogenity = 1-(goldEntropy-calMutual)/goldEntropy;
		
		System.out.println("v measure "+2*completeness*homogenity/(completeness+homogenity));
		
	//	System.out.println(1/((1/completeness)+(1/homogenity)));
		
		CloseFile(goldstandardFile);
		CloseFile(resultsFile);
	}
	
	public static double CalEntropy(Hashtable<String, Hashtable<String, Integer>> resultClusters){
		
		double entropy=0;
		
		Enumeration<String> k = resultClusters.keys();
		while (k.hasMoreElements()) {
			String cluster = (String) k.nextElement();
			double weight = (double)(resultClustersSize.get(cluster)/(double)corpusSize); //cluster weight
			
			double clusterEntropy = 0;
			Hashtable<String, Integer> clusterWords = resultClusters.get(cluster);
			Enumeration<String> words = clusterWords.keys();
			while (words.hasMoreElements()) {
				String nextWord = words.nextElement();
				int wordFreq = clusterWords.get(nextWord);
				double pxi = (double)wordFreq/(double)(resultClustersSize.get(cluster));
				//clusterEntropy += pxi *(Math.log(pxi)/Math.log(10));
				clusterEntropy += pxi *Math.log(pxi);
			}
			entropy += weight*clusterEntropy;
		}
		return (-1)*entropy;
	}
	
	public static double CalMutualInformation(){
		
		double mutualInf=0;
		
		Enumeration<String> results = resultClusters.keys();
		while (results.hasMoreElements()) {
			
			String c1 = results.nextElement();
			Enumeration<String> golds = goldClusters.keys();
			while (golds.hasMoreElements()) {
				
				String c2 = golds.nextElement();
				Hashtable<String, Integer> c1Words = resultClusters.get(c1);
				Hashtable<String, Integer> c2Words = goldClusters.get(c2);
				
				int intersection = Intersection(c1Words, c2Words);
				
				double pInter = (double)intersection/(double)counter;
				double pX = (double)resultClustersSize.get(c1)/(double)counter;
				double pY = (double)goldClustersSize.get(c2)/(double)counter;
				if(pInter!=0)
					//mutualInf += pInter*(Math.log(pInter/(pX*pY))/Math.log(10));
					mutualInf += pInter*Math.log(pInter/(pX*pY));
			}
		}
		
		return mutualInf;
	}
	
	public static int Intersection(Hashtable<String, Integer> c1Words, Hashtable<String, Integer> c2Words){
		
		int intersection=0;
		
		Enumeration<String> c2 = c2Words.keys();
		while (c2.hasMoreElements()) {
			
			String nextC2Word = c2.nextElement();
			if(c1Words.containsKey(nextC2Word)){
				int c2Int = c2Words.get(nextC2Word);
				int c1Int = c1Words.get(nextC2Word);
				if(c1Int<c2Int)
					intersection += c1Int;
				else
					intersection += c2Int;
			}
		}
		
		return intersection;
	}
	
	public static double CalEntropyGold(Hashtable<String, Hashtable<String, Integer>> goldClusters){
		
		double entropy=0;
		
		Enumeration<String> k = goldClusters.keys();
		while (k.hasMoreElements()) {
			String cluster = (String) k.nextElement();
			double weight = (double)(goldClustersSize.get(cluster)/(double)counter); //cluster weight
			
			double clusterEntropy = 0;
			Hashtable<String, Integer> clusterWords = goldClusters.get(cluster);
			Enumeration<String> words = clusterWords.keys();
			while (words.hasMoreElements()) {
				String nextWord = words.nextElement();
				int wordFreq = clusterWords.get(nextWord);
				double pxi = (double)wordFreq/(double)(goldClustersSize.get(cluster));
				clusterEntropy += pxi *Math.log(pxi);
			}
			entropy += weight*clusterEntropy;
		}
		return (-1)*entropy;
	}

	public static void AddToHashHash(String cluster, String word, Hashtable<String, Hashtable<String, Integer>> hashash){
		
		if(hashash.containsKey(cluster)){
			Hashtable<String, Integer> hash = hashash.get(cluster);
			AddToHashtable_str(hash,word);
		}
		else{
			Hashtable<String, Integer> newHash = new Hashtable<String, Integer>();
			AddToHashtable_str(newHash,word);
			hashash.put(cluster, newHash);
		}
	}
	
	public static void AddToHashtable_str(Hashtable<String, Integer> table, String key){
		
		if(table.containsKey(key)){
			int currentCount = table.get(key);
			table.put(key, currentCount+1);
		}
		else
			table.put(key, 1);
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
