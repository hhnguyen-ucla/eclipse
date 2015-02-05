package quant.checkquant;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

public class quantArea {
	//improvement of quantArea and quantArea2, for fractions (2), for frac correction(2),
	//for duplicated area (none)

	public static void main(String[] arg) throws IOException{
		TreeMap<String, TreeMap<String, Double>> infotree=new TreeMap<>();		
//		TreeMap<String, Double> load=new TreeMap<>();
//		load.put("2", (double)8/20);
//		load.put("3", (double)8/20); //total 20
//		load.put("4", (double)8/20);
//		load.put("5", (double)8/20);
//		load.put("6", (double)8/20);
//		load.put("7", (double)10/20);
//		load.put("8", (double)10/20);
//		load.put("9", (double)10/20);
//		load.put("1", (double)18/25); //total 25, 18/25*20
//		infotree.put("e:/wolfei/20141106-quantArea/ecoli/20141024ecoli.txt",load);
//		TreeMap<String, Double> load2=new TreeMap<>();
//		load2.put("2", (double)15/20);
//		load2.put("3", (double)15/20); //total 20
//		load2.put("4", (double)15/20);
//		load2.put("5", (double)15/20);
//		load2.put("6", (double)15/20);
//		load2.put("7", (double)15/20);
//		load2.put("8", (double)15/20);
//		load2.put("9", (double)15/20);
//		load2.put("1", (double)15/20); //total 25, 18/25*20
//		infotree.put("e:/wolfei/20141106-quantArea/ecoli/20150112ecoli.txt",load2);
		TreeMap<String, Double> load2=new TreeMap<>();
		load2.put("2", (double)8/25);
		load2.put("3", (double)8/25); //total 20
		load2.put("4", (double)8/25);
		load2.put("5", (double)8/25);
		load2.put("6", (double)8/25);
		load2.put("7", (double)10/25);
		load2.put("8", (double)10/25);
		load2.put("9", (double)10/25);
		load2.put("1", (double)14/27); //total 25, 18/25*20
		infotree.put("e:/wolfei/20141106-quantArea/quantArea2/20141024wolfeiC.txt",load2);
		infotree.put("e:/wolfei/20141106-quantArea/quantArea2/20141024wolfeiC_2.txt",load2);
		infotree.put("e:/wolfei/20141106-quantArea/quantArea2/20141024wolfeiB.txt",load2);
		TreeMap<String, Double> load3=new TreeMap<>();
		load3.put("2", (double)4/20);
		load3.put("3", (double)8/20); //total 20
		load3.put("4", (double)8/20);
		load3.put("5", (double)8/20);
		load3.put("6", (double)8/20);
		load3.put("7", (double)10/20);
		load3.put("8", (double)10/20);
		load3.put("9", (double)10/20);
		load3.put("1", (double)14/20); //total 25, 18/25*20
		infotree.put("e:/wolfei/20141106-quantArea/quantArea2/20141024wolfeimC.txt",load3);
		TreeMap<String, Double> load4=new TreeMap<>();
		load4.put("2", (double)4/20);
		load4.put("3", (double)4/20); //total 20
		load4.put("4", (double)5/20);
		load4.put("5", (double)5/20);
		load4.put("6", (double)5/20);
		load4.put("7", (double)8/20);
		load4.put("8", (double)8/20);
		load4.put("9", (double)8/20);
		load4.put("1", (double)10/20); //total 25, 18/25*20
		infotree.put("e:/wolfei/20141106-quantArea/quantArea2/20140928wolfeiB.txt",load4);
		infotree.put("e:/wolfei/20141106-quantArea/quantArea2/20140928wolfeiC.txt",load4);
		infotree.put("e:/wolfei/20141106-quantArea/quantArea2/20140930wolfeimC.txt",load4);
		
		quantArea test=new quantArea();
		for (Map.Entry<String, TreeMap<String, Double>> e:infotree.entrySet()){
			String path=e.getKey();
//			System.out.println(path);
			TreeMap<String, Double> tree=e.getValue();
			test.read(path,path.replace("txt", "out"), path.replace("txt", "temp"),true,tree);
		}
	}
	
	public void read(String path,String out,String temp,boolean all,TreeMap<String,Double>load) throws IOException{
		BufferedWriter bwtemp = new BufferedWriter(new FileWriter(temp));
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		BufferedReader br = new BufferedReader(new FileReader(path));
		String line=br.readLine();line=br.readLine();
//		System.out.println(path);
		TreeMap<String, ArrayList<Double>>temptree=new TreeMap<>();
		TreeMap<String, TreeMap<String, ArrayList<Double>>> tree=new TreeMap<>();
		while (line!=null){
			String[] array=line.split("\t");
			String protein=array[5];
			String peptide=array[1];
			String lastres=peptide.substring(peptide.length()-1);
//			System.out.println(lastres);
			int countk=0;
			for (int i=0;i<peptide.length();i++){
				String c=peptide.charAt(i)+"";
				if (c.equals("k"))
					countk++;
			}
			String mod=array[6].replace(";"	, ",");
			String conf=array[0];
			String area=array[13];
			double score=Double.parseDouble(array[14]);
			String sample=array[28].replace(".raw", "");
			String frac=whatfrac(sample);
			String key="";
			if (all==true)
				key=protein;
			else
				key=protein+"\t"+sample;
			String ambi=array[2];
			
			boolean check=false;
			if (conf.contains("High") && ambi.contains("Unambig")&& peptide.length()>7&& score>25&&!lastres.equals("k")&&countk<2){
				check=true;
			}
			
			if (check==true&& area.length()==0){
//				if (check==true&& area.length()==0 &&mod.contains("K")){
				String tkey=key+"\t"+peptide+"\t"+mod+"\t\t";
				if (temptree.containsKey(tkey)){
					ArrayList<Double>tvalue=temptree.get(tkey); //highest score - number of spec
					double tscore=tvalue.get(0);
					double tnum=tvalue.get(1)+1;
					if (tscore<score)
						tscore=score;
					ArrayList<Double> newvalue=new ArrayList<Double>();
					newvalue.add(tscore); newvalue.add(tnum);
					temptree.put(tkey, newvalue);
				}else{
					ArrayList<Double> newvalue=new ArrayList<Double>();
					newvalue.add(score); newvalue.add(1.0);
					temptree.put(tkey, newvalue);
				}
			}
			
			if (check==true && area.length()>0){	
				String ckey=peptide+";"+mod+":"+area+"\t"+frac; //remove peptides with the same areas
				if (tree.containsKey(key)){
					TreeMap<String, ArrayList<Double>> child=tree.get(key);
					if (child.containsKey(ckey)){
						ArrayList<Double> cvalue=child.get(ckey);
						double cscore=cvalue.get(0);
						double cnum=cvalue.get(1)+1;
						if (cscore<score)
							cscore=score;
						ArrayList<Double>newvalue=new ArrayList<Double>();
						newvalue.add(cscore); newvalue.add(cnum);
						child.put(ckey,newvalue);
					}else{
						ArrayList<Double>newvalue=new ArrayList<Double>();
						newvalue.add(score); newvalue.add(1.0);
						child.put(ckey, newvalue);
					}
					if (!key.equals(""))
						tree.put(key, child);
				}else{
					TreeMap<String, ArrayList<Double>> child=new TreeMap<>();					
					ArrayList<Double>newvalue=new ArrayList<Double>();
					newvalue.add(score); newvalue.add(1.0);
					child.put(ckey, newvalue);
					if (!key.equals(""))
						tree.put(key, child);
				}
			}
			line=br.readLine();
		}
		TreeMap<String,String> allprot=sumPep(tree, load); //reconcile peptides from different fractions
		//key=prot, value=pep;mod\tarea\tscore;numpep - 
		
		BufferedWriter b=new BufferedWriter(new FileWriter("e:/temp.txt"));
		//write prot vs pep vs area vs numspec
		for (Map.Entry<String, String> entry: allprot.entrySet()){ //loop each protein
			String[] array=entry.getValue().split("-");
//			System.out.println(entry.getValue());
			for (int i=0;i<array.length;i++){
				String[] a=array[i].split("\t");
				String donekey=entry.getKey()+"\t"+a[0].replace(";","\t")+"\t\t";
				String[] aa=a[2].split(";");
				double numspec=Double.parseDouble(aa[1]);
				if (temptree.containsKey(donekey)){
					ArrayList<Double>list=temptree.get(donekey); //[score,numpep]
					double tempnumspec=(list.get(1));
					numspec=numspec+tempnumspec;
					temptree.remove(donekey);
				}
				int spec=(int)Math.round(numspec);
				bwtemp.write(entry.getKey()+"\t"+array[i].replace(";","\t").replace(aa[1],spec+"")+"\n");
			}
		}
		
		for (Map.Entry<String, ArrayList<Double>> t:temptree.entrySet()){
//			b.write("thisis "+t.getKey()+" temptree\n");
			ArrayList<Double>l=t.getValue();
			int numspec=(int)Math.round(l.get(1));
			bwtemp.write(t.getKey()+l.get(0)+"\t"+numspec+"\n");
		}
		
		TreeMap<String, Double> quanttree=quant(allprot);
		System.out.println(tree.size()+" "+quanttree.size());
		for (Map.Entry<String, Double> entry:quanttree.entrySet()){
			double area=entry.getValue();
			double log=Math.log10(area);
			bw.write(entry.getKey()+"\t"+entry.getValue()+"\t"+log+"\n");
		}
		br.close();bw.flush();bw.close(); bwtemp.flush();bwtemp.close();
		b.flush();b.close();
	}
	
	private TreeMap<String, Double> quant(TreeMap<String, String> tree) throws IOException{
		TreeMap<String, Double> returns=new TreeMap<>();
		int count=0;
		for (Map.Entry<String, String> entry: tree.entrySet()){ //loop each protein
			String[] array=entry.getValue().split("-");

			
			//calculate ave area of top 3, if not have 3: just the rest
			int top=5;
			if (array.length>=top){
				TreeMap<Long,String> peps=new TreeMap<>();
				ArrayList<Long>list=new ArrayList<>(); //inten list

				for (int i=0;i<array.length;i++){ //loop each peptide;mod
//					System.out.println("from quant "+array[i]);
					String[] a=array[i].split("\t"); //get inten for each peptide:mod
						long inten=(long)Math.round(Double.parseDouble(a[1]));
						list.add(inten);
						peps.put(inten,a[0]);
				}
				Collections.sort(list);
				int n=list.size();
				ArrayList<Long> inten=new ArrayList<>(); double topn=0; String pepn="";
				for (int i=1;i<=top;i++){
					long nth=(list.get(n-i));
					topn+=list.get(n-i);
					pepn+=peps.get(nth)+":"+nth+" ,";
				}
				returns.put(entry.getKey()+"\t"+array.length+"\t"+pepn, (double)Math.round(topn/n));
			}else{
				double topn=0;
				String pepn=entry.getValue().replace("-", ",").replace("\t", ":");
				for (int i=0;i<array.length;i++){ //loop each peptide;mod
					String[] a=array[i].split("\t"); //get inten for each peptide:mod
					double inten=Double.parseDouble(a[1]);
					topn+=inten;
					
				}
				returns.put(entry.getKey()+"\t"+array.length+"\t"+pepn, (double)Math.round(topn/array.length));
			}
		}		
		return returns;
	}
	
	//work across fractions
	private TreeMap<String, String> sumPep(TreeMap<String, TreeMap<String,ArrayList<Double>>> tree,TreeMap<String, Double>load){
		TreeMap<String, String> returns=new TreeMap<>();
		for (Map.Entry<String, TreeMap<String,ArrayList<Double>>> entry: tree.entrySet()){ //loop each protein
			String key=entry.getKey(); //prot
			TreeMap<String,ArrayList<Double>> peps=entry.getValue();//ckey=peptide+";"+mod+":"+area+"\t"+frac;value=[score,numpep]
			//key=pep;mod,area\tfrac value=score, numpep
			TreeMap<String, ArrayList<Double>> allpep=new TreeMap<>(); //pep;mod vs area,score,numpep
			for (Map.Entry<String, ArrayList<Double>> pscore:peps.entrySet()){
				String pep=pscore.getKey();
				ArrayList<Double>pvalue=pscore.getValue(); //score, num spec
				double score=pvalue.get(0);
				double numpep=pvalue.get(1);
				String[] array=pep.replace(":", "\t").split("\t");
				double area=Double.parseDouble(array[1]);
				double ratio=load.get(array[2]);
				area=area/ratio;
				if (allpep.containsKey(array[0])){
					ArrayList<Double>list=allpep.get(array[0]);
					double newarea=area+list.get(0);
					double newscore=list.get(1);
					double newnum=list.get(2)+numpep;
					if (newscore<score)
						newscore=score;
					ArrayList<Double>newlist=new ArrayList<>();
					newlist.add(newarea); newlist.add(newscore);newlist.add(newnum);
					allpep.put(array[0],newlist);
				}else{
					ArrayList<Double> newlist=new ArrayList<>();
					newlist.add(area); newlist.add(score);newlist.add(numpep);
					allpep.put(array[0], newlist);
				}
			}
			String value=""; //pep;mod\tarea\tscore;numpep - 
			for (Map.Entry<String, ArrayList<Double>> e:allpep.entrySet()){
				ArrayList<Double>list=e.getValue();
				value+=e.getKey()+"\t"+list.get(0)+"\t"+list.get(1)+";"+list.get(2)+"-";
			}
			returns.put(key, value);
		}
		return returns;
	}
	
	private String  whatfrac(String sample){
		String returns="";
		if (sample.contains("2."))
			returns="2";
		else if (sample.contains("3."))
			returns="3";
		else if (sample.contains("4."))
			returns="4";
		else if (sample.contains("5."))
			returns="5";
		else if (sample.contains("6."))
			returns="6";
		else if (sample.contains("7."))
			returns="7";
		else if (sample.contains("8."))
			returns="8";
		else if (sample.contains("9."))
			returns="9";
		else
			returns="1";
		return returns;
	}	
}
