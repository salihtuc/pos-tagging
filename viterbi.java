
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;

public class viterbi {

	private static HashMap<String, NodeV> transitions;
	private static HashMap<String, NodeV> emissions;
	private static ArrayList<String> statess;
	private static HashSet<String> wordsAll;

	private static double[][] viterbi;
	private static int[][] vitPrev;
	static int states = 0;
	static int smoothed = 0;

	public static void main(String[] args) throws IOException {

		transitions = new HashMap<String, NodeV>();
		emissions = new HashMap<String, NodeV>();
		statess = new ArrayList<String>();
		wordsAll = new HashSet<String>();
		long start = 0;
		long elapsedTime = 0;

		System.out.println("=== Building The Model");
		start = System.nanoTime();
		buildModel();
		elapsedTime = System.nanoTime() - start;
		System.out.println("=== Building The Model Completed");
		System.out.println("Total time: " + elapsedTime);
		System.out.println();

		states = statess.size();

		System.out.println("=== Tagging the text");
		start = System.nanoTime();
		tagging("PennCorpus12.txt","denemeOut.txt");
		elapsedTime = System.nanoTime() - start;
		System.out.println("=== Tagging Completed");
		System.out.println("Total time: " + elapsedTime);
	}

	private static void leplaceSmoothing(String word) {
		Iterator it = emissions.entrySet().iterator();
		while (it.hasNext()) {
			Map.Entry pair = (Map.Entry) it.next();
			NodeV tmp = (NodeV) pair.getValue();

			HashMap<String, Double> tmp2 = tmp.vals;
			Iterator it2 = tmp2.entrySet().iterator();
			int countwr = tmp2.size();
			double newcount = tmp.val + countwr + 1;

			while (it2.hasNext()) {
				Map.Entry pair2 = (Map.Entry) it2.next();
				tmp2.replace((String) pair2.getKey(), ((double) pair2.getValue() * tmp.val) / newcount);
			}
			tmp2.put(word, (1.0 / newcount));
			wordsAll.add(word);
			tmp.val = newcount;
		}
		// System.out.println("\"" + word + "\" has been smoothed");
		smoothed++;
	}

	private static void tagging(String sysinput,String sysoutput) throws IOException {
		int countWord = 0;
		int correctTag = 0;
		int cs = 0;

		BufferedWriter output = null;
		BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(sysinput), "UTF8"));

		try {

			File file = new File(sysoutput);
			output = new BufferedWriter(new FileWriter(file));

			String sentence = br.readLine();

			while (sentence != null) {
				sentence="<s> "+sentence;
						String[] words=sentence.toLowerCase().split("\\s+");
						viterbi = new double[states][words.length];
						vitPrev = new int[states][words.length];

						for (int i = 0; i < words.length; i++) {

							String word = words[i];
							if (!wordsAll.contains(word)) {
//								System.out.println(word+" smoothed");
								leplaceSmoothing(word);
							}

							for (int j = 0; j < statess.size(); j++) {
//								System.out.println(statess.get(j)+" "+word);
								double score = vitScore(statess.get(j), word, i, j);
//								System.out.println(statess.get(j)+" "+word+" "+score);
								viterbi[j][i] = score;
							}

						}

						String tag = "";
						int idx = 0;
						double max = 0;
						int n = words.length - 1;
						for (int i = 0; i < states; i++) {
							if (viterbi[i][n] > max) {
								max = viterbi[i][n];
								idx = i;
							}
						}
						for (int i = n; i >= 1; i--) {

							String cek = statess.get(idx);
//							String truetag = getTag(words[i]);
							countWord++;
//							if (cek.equals(truetag)) {
//								correctTag++;
//							}

							tag = words[i] + "/" + cek + " " + tag;
							
							idx = vitPrev[idx][i];

						}
						output.write(tag + "\n");
//						words.clear();
						cs++;	
					

					sentence = br.readLine();
				
			}
		} finally {
			br.close();
			if (output != null) {
				output.close();
			}
		}

		System.out.println("Total evaluate sentence: " + cs);
		System.out.println("Total smoothed words: " + smoothed);
		System.out.println("Total evaluate words: " + countWord);
//		System.out.println("Total correct tag: " + correctTag);
//		System.out.println("Correctness persentage: " + (((double) correctTag / countWord) * 100) + "%");

	}

	private static double vitScore(String tag, String word, int idxWord, int idxTag) {
		double max = 0.0;
		double emProb = 0.0;
		NodeV emtag = emissions.get(tag);
		if (emtag != null) {
			if (emtag.vals.containsKey(word)) {
				emProb = emtag.vals.get(word);
//				System.out.println(word+" "+tag+" "+emtag.vals.get(word));
			}
		}

		if (idxWord == 0) {
			double transProb = 0.0;
			NodeV trans = transitions.get("<s>");
			if (trans.vals.containsKey(tag)) {
				transProb = trans.vals.get(tag);
			}
			max = (1 * transProb * emProb);
			vitPrev[idxTag][idxWord] = 0;
		} else {

			for (int i = 0; i < statess.size(); i++) {
				String prevTag = statess.get(i);
				double prevVitS = viterbi[i][idxWord - 1];
				double transProb = 0.0;

				NodeV trans = transitions.get(prevTag);
				if (trans != null) {
					if (trans.vals.containsKey(tag)) {
						transProb = trans.vals.get(tag);
					}
				}

				double score = (prevVitS * transProb * emProb);
				if (score > max) {
					max = score;
					vitPrev[idxTag][idxWord] = i;
				}
			}
		}
		return max;
	}

	private static void buildModel() throws IOException {

		BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream("out.txt"), "UTF8"));
		String tagPrev = "";
		String tag = "";
		String word = "";
		try {
			String line = br.readLine();
			int i=1;
			while (line != null) {
				if (i<170){
					String[] words = line.split(" ");
					tag = words[0].split("\\|")[1];
					tagPrev = words[0].split("\\|")[0];
					addTag(tagPrev, tag,Double.parseDouble(words[1]));
						
						
				}
				else{
					String[] words = line.split(" ");
					tag = words[0].split("\\|")[0];
					word = words[0].split("\\|")[1];
					addWord(tag, word,Double.parseDouble(words[1]));
					wordsAll.add(word);						
					}
					
				i++;
				line = br.readLine();
			}

		} finally {
			br.close();
		}
//		System.out.println(wordsAll.size());
		
//		for(String a:transitions.keySet()){
//			System.out.println(a+" "+transitions.get(a).name);
//			HashMap<String, Double> vals=transitions.get(a).vals;
//			for(String b:vals.keySet()){
//				System.out.println(b+" "+vals.get(b));
//			}
//		}

	}

	private static void addTag(String prev, String tag,Double value) {
		if (transitions.containsKey(prev)) {
			NodeV check = transitions.get(prev);
			check.addVal(tag,value);
		} else {
			NodeV check = new NodeV(prev, tag,value);
			transitions.put(prev, check);
			statess.add(prev);
		}
	}

	private static void addWord(String tag, String word,Double value) {
//		word = word.replaceAll("\t", "");
		if (!wordsAll.contains(word)) {
			wordsAll.add(word);
		}
		if (emissions.containsKey(tag)) {
			NodeV check = emissions.get(tag);
			check.addVal(word,value);
		} else {
			NodeV check = new NodeV(tag, word,value);
			emissions.put(tag, check);
		}
	}

	private static String getTag(String str) {
		int idx = str.lastIndexOf('/') + 1;
		return str.substring(idx);
	}

	private static String getWord(String str) {
		int idx = str.lastIndexOf('/');
		return str.substring(0, idx);
	}
}

class NodeV {

	String name;
	HashMap<String, Double> vals;
	double val;

	public NodeV(String name, String firstVal,double vall) {
		this.name = name;
		this.val = vall;
		vals = new HashMap<String, Double>();
		addVal(firstVal,vall);
	}

	public void addVal(String val,double vall) {
//		System.out.println(val+" "+vall);
//		if (vals.containsKey(val)) {
//			double curr = vals.get(val) +vall;
//			vals.replace(val, curr);
//		} else {
			vals.put(val, vall);
//			System.out.println(val+" "+vall);
//		}
		
	}

}
