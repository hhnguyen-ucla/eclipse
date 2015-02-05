package quant;

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

public class IPquant {
//improve fractions: put sample name at the end
	//improve temp: if modified peptide has both area and 0, delete 0

	public static void main(String[] arg) throws IOException{
		IPquant test=new IPquant();
		String path="e:/wolfei/20150113-modPeps/quantArea/IP-wolfeiBAcfracs.txt";
		test.read(path,path.replace("txt", "out"), path.replace("txt", "temp"),true,false); //ID=protein
		String path2="e:/wolfei/20150113-modPeps/quantArea/IP-wolfei.txt";
		test.read(path2,path2.replace("txt", "out"), path2.replace("txt", "temp"),false,false); //ID=protein+sample
	}
	
	public void read(String path,String out,String temp,boolean all,boolean confi) throws IOException{
		BufferedWriter bwtemp = new BufferedWriter(new FileWriter(temp));
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		BufferedReader br = new BufferedReader(new FileReader(path));
		String line=br.readLine();line=br.readLine();
//		System.out.println(line);
		TreeMap<String, Double>temptree=new TreeMap<>();
		TreeMap<String, TreeMap<String, Double>> tree=new TreeMap<>();
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
			if (ambi.contains("Unambig")&& peptide.length()>7&& score>25&&!lastres.equals("k")&&countk<2){
				if (confi==true){
					if (conf.contains("High"))
						check=true;
				}else
					check=true;
			}
			
			if (check==true&& area.length()==0 &&mod.contains("K")){
				String tkey=key+"\t"+peptide+"\t"+mod+"\t\t";
				if (temptree.containsKey(tkey)){
					double tscore=temptree.get(tkey);
					if (tscore<score)
						tscore=score;
					temptree.put(tkey, tscore);
				}else
					temptree.put(tkey, score);
			}
			
			if (check==true && area.length()>0){	
				String ckey=peptide+";"+mod+":"+area+"\t"+frac;
				if (tree.containsKey(key)){
					TreeMap<String, Double> child=tree.get(key);
					if (child.containsKey(ckey)){
						double cscore=child.get(ckey);
						if (cscore<score)
							child.put(ckey, score);
					}else{
						child.put(ckey, score);
					}
					if (!key.equals(""))
						tree.put(key, child);
				}else{
					TreeMap<String, Double> child=new TreeMap<>();
					child.put(ckey, score);
					if (!key.equals(""))
						tree.put(key, child);
				}
			}
			line=br.readLine();
		}
		TreeMap<String,String> allprot=sumPep(tree); //remove peps with the same areas
		
		//write prot vs pep vs area
		for (Map.Entry<String, String> entry: allprot.entrySet()){ //loop each protein
				String[] array=entry.getValue().split("-");
				for (int i=0;i<array.length;i++){
					if (all==false){
						String[] a=entry.getKey().split("\t");
						String sample=a[1];
						bwtemp.write(entry.getKey().replace("\t"+sample, "")+"\t"+array[i].replace(";","\t")+"\t"+sample+"\n");				
					}else					
						bwtemp.write(entry.getKey()+"\t"+array[i].replace(";","\t")+"\n");
				}
			if (temptree.containsKey(entry.getKey())){
				temptree.remove(entry.getKey()); //delete peptide entry if it has an area already
				System.out.println(entry.getKey());
			}
		}
		
		for (Map.Entry<String, Double> t:temptree.entrySet()){
			if (all==false){
				String[] a=t.getKey().split("\t");
				String sample=a[1];
				bwtemp.write(t.getKey().replace("\t"+sample, "")+t.getValue()+"\t"+sample+"\n");				
			}else
				bwtemp.write(t.getKey()+t.getValue()+"\n");
		}
		
		TreeMap<String, Double> quanttree=quant(allprot);
		System.out.println(tree.size()+" "+quanttree.size());
		for (Map.Entry<String, Double> entry:quanttree.entrySet()){
			double area=entry.getValue();
			double log=Math.log10(area);
			bw.write(entry.getKey()+"\t"+entry.getValue()+"\t"+log+"\n");
		}
		br.close();bw.flush();bw.close(); bwtemp.flush();bwtemp.close();
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
	
	private TreeMap<String, String> sumPep(TreeMap<String, TreeMap<String,Double>> tree){
		TreeMap<String, String> returns=new TreeMap<>();
		for (Map.Entry<String, TreeMap<String,Double>> entry: tree.entrySet()){ //loop each protein
			String key=entry.getKey();
			TreeMap<String,Double> peps=entry.getValue(); //pep;mod:area\tfrac
			TreeMap<String, ArrayList<Double>> allpep=new TreeMap<>(); //pep;mod vs area,score
			for (Map.Entry<String, Double> pscore:peps.entrySet()){
				String pep=pscore.getKey();
				double score=pscore.getValue();
				String[] array=pep.replace(":", "\t").split("\t");
				double area=Double.parseDouble(array[1]);
				if (allpep.containsKey(array[0])){
					ArrayList<Double>list=allpep.get(array[0]);
					double newarea=area+list.get(0);
					double newscore=list.get(1);
					if (newscore<score)
						newscore=score;
					ArrayList<Double>newlist=new ArrayList<>();
					newlist.add(newarea); newlist.add(newscore);
					allpep.put(array[0],newlist);
				}else{
					ArrayList<Double> newlist=new ArrayList<>();
					newlist.add(area); newlist.add(score);
					allpep.put(array[0], newlist);
				}
			}
			String value=""; //pep;mod\tarea - 
			for (Map.Entry<String, ArrayList<Double>> e:allpep.entrySet()){
				ArrayList<Double>list=e.getValue();
				value+=e.getKey()+"\t"+list.get(0)+"\t"+list.get(1)+"-";
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
