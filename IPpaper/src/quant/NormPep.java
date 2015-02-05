package quant;

import java.io.BufferedReader;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.TreeMap;

public class NormPep {

	public static void main(String[] arg) throws IOException{
		NormPep test=new NormPep();
//		test.process();
		String path="e:/wolfei/20141106-quantArea/quantArea2/histogramModPep.txt";
		BufferedWriter bw = new BufferedWriter(new FileWriter(path));
		String[] paths={"e:/wolfei/20141106-quantArea/quantArea2/20140928wolfeiB.temp2","e:/wolfei/20141106-quantArea/quantArea2/20140928wolfeiC.temp2","e:/wolfei/20141106-quantArea/quantArea2/20140930wolfeimC.temp2","e:/wolfei/20141106-quantArea/quantArea2/20141024wolfeiB.temp2","e:/wolfei/20141106-quantArea/quantArea2/20141024wolfeiC.temp2","e:/wolfei/20141106-quantArea/quantArea2/20141024wolfeiC_2.temp2","e:/wolfei/20141106-quantArea/quantArea2/20141024wolfeimC.temp2"};
		String[] modlist={"Acetyl","K+70","K+86","K+84","K+68"};	
		for (int i=0;i<paths.length;i++){
			for (int j=0;j<modlist.length;j++){
				bw.write(modlist[j]+"\n"+test.histograms(paths[i],modlist[j])+"\n");
			}
		}
		bw.flush();bw.close();
	}
	
	public String histograms(String path,String mod) throws IOException{
		String returns=path+"\n";
		BufferedReader br = new BufferedReader(new FileReader(path));
		String line=br.readLine();
		ArrayList<Double> loglist=new ArrayList<Double>();
		while (line!=null){
			String[] array=line.split("\t");
			String areatxt=array[3];
			if (line.contains(mod)&&!areatxt.equals("")){
				double area=Double.parseDouble(areatxt);
				loglist.add(Math.log10(area));
			}
			line=br.readLine();
		}
		String[] range={"6.0-6.5","6.5-7.0","7.0-7.5","7.5-8.0","8.0-8.5","8.5-9.0","9.0-9.5","9.5-10.0","10.0-10.5","10.5-11.0","11.0-11.5","11.5-12","12.0-12.5","12.5-13.0"};
		TreeMap<String, Integer> histotree=makeHistogram(loglist);
		for (int i=0;i<range.length;i++){
			String count=0+"";
			if (histotree.containsKey(range[i]))
				count=histotree.get(range[i])+"";
			returns+=count+"\n";
//			returns+=range[i]+"\t"+count+"\n";
		}
		return returns;
	}
	
	private TreeMap<String,Integer> makeHistogram(ArrayList<Double>list){	
		TreeMap<String,Integer> count=new TreeMap<String, Integer>();
		for (int i=0;i<list.size();i++){
			double area=list.get(i);
			double start=6.0;
			while (start+0.5<13.0){
				if (area>=start && area<start+0.5){
					String range=start+"-"+(double)(start+0.5);
					if (count.containsKey(range)){
						count.put(range, count.get(range)+1);
					}else
						count.put(range,1);
				}
				start=start+0.5;
			}
		}
		return count;
	}
	
	//normalize peptide area based on total protein area (sum of ave top 5 as in AveProt)
	//work with temp, write temp2
	public void normpep(String path,String out,double ratio) throws IOException{
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		BufferedReader br = new BufferedReader(new FileReader(path));
		String line=br.readLine();
		double sum=0;
		while (line!=null){
			String[]array=line.split("\t");
			String areatxt=array[3];
			String newarea="";
			if (!areatxt.equals("")){
				double area=Double.parseDouble(areatxt);
				area=area*ratio;
				newarea=area+"";
				sum+=area;
			}
			bw.write(line.replace(areatxt, newarea)+"\n");	
			line=br.readLine();
		}
		System.out.println(sum+"\t"+path);
		br.close();bw.flush();bw.close();
	}
	void process() throws IOException{
		String allprotout="e:/wolfei/20141106-quantArea/allwolfei.out";
		int samplenum=7; double maxarea=0;
		ArrayList<Double> sumlist=new ArrayList<>();
		for (int i=0;i<samplenum;i++){
			double sum=calSum(allprotout,i*4+3,i*4+1,0); //String path,int areacol,int pepcol,int protcol
			if (sum>maxarea)
				maxarea=sum;
			sumlist.add(sum);
		}
		
		String[] paths={"e:/wolfei/20141106-quantArea/quantArea2/20140928wolfeiB.temp","e:/wolfei/20141106-quantArea/quantArea2/20140928wolfeiC.temp","e:/wolfei/20141106-quantArea/quantArea2/20140930wolfeimC.temp","e:/wolfei/20141106-quantArea/quantArea2/20141024wolfeiB.temp","e:/wolfei/20141106-quantArea/quantArea2/20141024wolfeiC.temp","e:/wolfei/20141106-quantArea/quantArea2/20141024wolfeiC_2.temp","e:/wolfei/20141106-quantArea/quantArea2/20141024wolfeimC.temp"};
		for (int i=0;i<7;i++){
			String path=paths[i];
			double sum=sumlist.get(i);
			double ratio=maxarea/sum;
			normpep(path,path+"2",ratio);
//			System.out.println(ratio+" "+path);
		}
	}
	
	//at least 2 pep
	private double calSum(String path,int areacol,int pepcol,int protcol) throws IOException{
		double sum=0;
		BufferedReader br = new BufferedReader(new FileReader(path));
		String line=br.readLine();
		while (line!=null){
			String[] array=line.replace("\t\t", "\t-\t-").split("\t");
			String pep=array[pepcol].replace("-", "");
			String area=array[areacol].replace("-", "");
			String prot=array[protcol];
			if (!prot.contains("P")&&!prot.contains("Q")&&!prot.contains(";")){
				if (pep.length()>0){
					int pepnum=Integer.parseInt(pep);
					if (pepnum>2){
						if (area.length()>0){
							double areanum=Double.parseDouble(area);
							sum+=areanum;
						}
					}
				}
			}
			line=br.readLine();
		}
		br.close();
		return sum;
	}

}

