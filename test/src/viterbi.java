
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
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

	private static HashMap<String, Double> transitions;
	private static HashMap<String, Double> emissions;
	private static HashMap<String, Double> tagValues;
	private static HashMap<String, String> tagConvert;
	private static HashMap<String, Double> tagValuesEmission;
	private static ArrayList<String> statess;
	private static HashSet<String> wordsAll;

	private static double[][] viterbi;
	private static int[][] vitPrev;
	static int states = 0;
	static int smoothed = 0;

	public static void main(String[] args) throws IOException {

		transitions = new HashMap<String, Double>();
		emissions = new HashMap<String, Double>();
		tagConvert = new HashMap<String, String>();
		
		tagValues = new HashMap<String, Double>();
		tagValuesEmission = new HashMap<String, Double>();
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
		// System.out.println(statess);
		
		System.out.println("=== Tagging the text");
		start = System.nanoTime();

		tagging("pennUntagged500.txt", "taggedp.txt");	// (UntaggedFile, the file which will tagged)
		elapsedTime = System.nanoTime() - start;
		System.out.println("=== Tagging Completed");
		System.out.println("Total time: " + elapsedTime);
	}

	private static void leplaceSmoothing(String word, String i) {
		Iterator it = emissions.entrySet().iterator();
		double sum = 0;
		while (it.hasNext()) {
			Map.Entry pair = (Map.Entry) it.next();
			String wordd = (String) pair.getKey();
			String tag = wordd.split("\\|", 2)[0];
			double value = (double) pair.getValue();
			if (tag.equals(i)) {
				sum += value + Double.MIN_VALUE;
//				sum += value + -1*Double.MAX_VALUE;
			}
		}

		while (it.hasNext()) {
			Map.Entry pair = (Map.Entry) it.next();
			String wordd = (String) pair.getKey();
			String tag = wordd.split("\\|", 2)[0];
			if (tag.equals(i)) {
				emissions.replace(wordd, emissions.get(wordd) / sum);
				//// sum += value + 0.00001;
			}
		}

		emissions.put(i + "|" + word, Double.MIN_VALUE / sum);
//		emissions.put(i + "|" + word, -1*Double.MAX_VALUE / sum);
		wordsAll.add(word);

		System.out.println("\"" + word + "\" has been smoothed");
		smoothed++;
	}

	private static void tagging(String sysinput, String sysoutput) throws IOException {
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
				sentence = sentence +" <end>";
				String[] words = sentence.toLowerCase().split("\\s+");
				viterbi = new double[states][words.length];
				vitPrev = new int[states][words.length];

				for (int i = 0; i < words.length; i++) {

					String word = words[i];
					if (!wordsAll.contains(word)) {
						//// System.out.println(word+" smoothed");
						leplaceSmoothing(word, String.valueOf(i));
					}

					for (int j = 0; j < statess.size(); j++) {
						//// System.out.println(statess.get(j)+" "+word);
						double score = vitScore(statess.get(j), word, i, j);
						//// System.out.println(statess.get(j)+" "+word+"
						//// "+score);
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
//				tag = "<end>" + "/" + statess.get(idx) + " " + tag;
				
				idx = vitPrev[idx][n];
				for (int i = n-1; i >= 0; i--) {

					String cek = statess.get(idx);
					//// String truetag = getTag(words[i]);
					countWord++;
					//// if (cek.equals(truetag)) {
					//// correctTag++;
					//// }

					tag = words[i] + "/" + cek + " " + tag;

					idx = vitPrev[idx][i];

				}
				output.write(tag + "\n");
				//// words.clear();
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
		//// System.out.println("Total correct tag: " + correctTag);
		//// System.out.println("Correctness persentage: " + (((double)
		//// correctTag /countWord) * 100) + "%");

	}

	private static double vitScore(String tag, String word, int idxWord, int idxTag) {
		double max = 0.0;
		double emProb = Double.MIN_VALUE;
//		double emProb = -1*Double.MAX_VALUE;
		String emit = tag + "|" + word;
		if (emissions.keySet().contains(emit))
			emProb = emissions.get(emit);
		
			

		if (idxWord == 0) {
			double transProb = Double.MIN_VALUE;
//			double transProb = -1*Double.MAX_VALUE;
			String trans = "<s>|" + tag;
			if (transitions.keySet().contains(trans))
				transProb = transitions.get(trans);
			//// NodeV trans = transitions.get("<s>");
			//// if (trans.vals.containsKey(tag)) {
			//// transProb = trans.vals.get(tag);
			//// }
			max = (1 * transProb * emProb);
			vitPrev[idxTag][idxWord] = 0;
		} else {

			for (int i = 0; i < statess.size(); i++) {
				String prevTag = statess.get(i);
				double prevVitS = viterbi[i][idxWord - 1];
				double transProb = Double.MIN_VALUE;
//				double transProb = -1*Double.MAX_VALUE;
				String trans = prevTag + "|" + tag;
				if (transitions.keySet().contains(trans))
					transProb = transitions.get(trans);
				
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

		BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream("outp.txt"), "UTF8"));
		String tagPrev = "<s>";
		String tag = "";
		String word = "";
		try {
			String line = br.readLine();
			int i = 1;
			while (line != null) {
				if (line.length() > 0) {
//					System.out.println(line);
					if (i <= 323) {
						String[] words = line.split(" ");
						tag = words[0].split("\\|", 2)[1];
						transitions.put(words[0], Double.parseDouble(words[1]));
//						transitions.put(words[0], -1 *Double.parseDouble(words[1]));
						if (!statess.contains(tag)) {
							statess.add(tag);
						}

					} else {
						String[] words = line.split(" ");
						word = words[0].split("\\|", 2)[1];
						
						wordsAll.add(word);
						emissions.put(words[0], Double.parseDouble(words[1]));
//						emissions.put(words[0], -1 *Double.parseDouble(words[1]));
					}

					i++;
					line = br.readLine();
				}

			}
			
			wordsAll.add("<start>");
			wordsAll.add("<end>");
			emissions.put("<s>|<start>", 1.0);
			emissions.put("</s>|<end>", 1.0);
//			emissions.put("<s>|<start>", -1.0);
//			emissions.put("</s>|<end>", -1.0);
			
		} finally {
			br.close();
		}
//		System.out.println(statess);

	}

	private static void calculate(String word, String tag, String prevTag) {
		
		String keyTrans = prevTag + "|" + tag;
		String keyEmit = tag + "|" + word;
		
		
		if (transitions.keySet().contains(keyTrans)) {
			double value = transitions.get(keyTrans);
			transitions.put(keyTrans, value + 1);
			tagValues.put(prevTag, tagValues.get(prevTag) + 1);
		} else {
			transitions.put(keyTrans, 1.0);
			tagValues.put(prevTag, 1.0);
		}
//		if (tagValuesEmission.keySet().contains(tag)) {
//			double value = tagValuesEmission.get(tag);
//			tagValuesEmission.put(tag, tagValuesEmission.get(tag) + 1);
//			
//		} else {
//			tagValuesEmission.put(tag, 1.0);
//		}
		if (!word.equals("<end>")){
		if (emissions.keySet().contains(keyEmit)) {
			double value = emissions.get(keyEmit);
			emissions.put(keyEmit, value + 1);
			
		} else {
			emissions.put(keyEmit, 1.0);
		}
		}
		
	}

	private static void avg() throws FileNotFoundException {
//		System.out.println(tagValues.get("</s>"));
		PrintWriter pw= new PrintWriter(new File("outp.txt"));	// Probabilities
		for (String key : transitions.keySet()) {
			String prevTag = key.split("\\|")[0];
//			System.out.println(prevTag);
			transitions.put(key, transitions.get(key) / tagValues.get(prevTag));
			pw.println(key+"_t "+transitions.get(key) / tagValues.get(prevTag));
		}

		for (String key : emissions.keySet()) {
			String tag = key.split("\\|")[0];
//			System.out.println(key);
			emissions.put(key, emissions.get(key) / tagValues.get(tag));
			pw.println(key+" "+emissions.get(key) / tagValues.get(tag));
		}
		
		emissions.put("<s>|<start>",1.0);
		emissions.put("</s>|<end>",1.0);
		
		pw.println("<s>|<start>"+" "+1.0);
		pw.println("</s>|<end>"+" "+1.0);
		
		pw.close();

	}

	private static void buildModel1() throws IOException {

		BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream("17tagged_500"), "UTF8"));
		
		String tag = "";
		String word = "";
		int tagCounter=1;
		try {
			String line = br.readLine();
			int i = 1;
			while (line != null) {
				if (line.length() > 0) {
					String tagPrev = "<s>";
					String[] words = line.split("\\s+");
					// System.out.println(sentence+" "+words.length+"
					// "+words[1]+" "+words[4]);
					for (String wordd : words) {
						int index = wordd.lastIndexOf("/");
//						tag = "t" + (Integer.parseInt(wordd.substring(index+1))+1);
						tag = "t" + wordd.substring(index+1);
//						System.out.println(tag);
//						if (tagConvert.keySet().contains(tag.trim())){
//							tag=tagConvert.get(tag.trim());
//						}
//						else{
//							tagConvert.put(tag.trim(), "t"+tagCounter);
//							tag="t"+tagCounter;
//							tagCounter++;
//						}
						word = wordd.substring(0, index).toLowerCase();
						
//						System.out.println(tag+"-"word);
						calculate(word, tag, tagPrev);
						tagPrev = tag;
					}
					calculate("<end>", "</s>", tagPrev);
					
					
				} 
//				tagValues.put("</s>", 1.0);
				if (tagValues.keySet().contains("</s>")) {
//					double value = tagValues.get(tag);
					tagValues.put("</s>", tagValues.get("</s>") + 1);

				} else {
					tagValues.put("</s>", 1.0);
				}
				line = br.readLine();
			}

		} finally {
			br.close();
		}
		avg();
		System.out.println(tagConvert.keySet());
		System.out.println(tagConvert.values());
	}


}
