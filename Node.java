import java.util.ArrayList;
import java.util.HashMap;

public class Node {
	int stateNum;
	String word;
	
	ArrayList<Integer> next;
	ArrayList<Integer> prev;
	
	ArrayList<Double> alpha;
	ArrayList<Double> beta;

	//HashMap<String, Double> tagScores;
	ArrayList<Double> tagScores;
	
	boolean isEndState;
	boolean isStartState;
	
	public Node(){
		
	}
	
	public Node(int stateNum, String word) {
		this.stateNum = stateNum;
		this.word = word;
		
		this.next = new ArrayList<>();
		this.prev = new ArrayList<>();
		
		this.alpha = new ArrayList<>();
		this.beta = new ArrayList<>();
		
		//this.tagScores = new HashMap<>();
		this.tagScores = new ArrayList<>();
		
		if(stateNum == 0){
			this.isStartState = true;
		}
		else{
			this.isStartState = false;
		}
		
		this.isEndState = false;
	}

	@Override
	public String toString() {
		return stateNum + ": " + word + " (next: " + next + " | prev: " + prev + " | is end state: " + isEndState + ") ";
	}
	
}
