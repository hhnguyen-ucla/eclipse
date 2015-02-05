package quant;

import java.io.BufferedReader;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

import jdk.jfr.events.FileWriteEvent;

public class countModPep {

	public static void main(String[] arg) throws IOException{
		countModPep test=new countModPep();
//		test.process1();
		test.process2();
	}
	void process2() throws IOException{
		File folder=new File("E:/wolfei/20141106/quantArea");
		File[] list=folder.listFiles();
		String[] paths=new String[list.length];
		for (int i=0;i<list.length;i++){
			String file=list[i].getAbsolutePath();
			paths[i]=file;
		}
		
		String out="e:/wolfei/20150113-modPeps/ReplicatedPep.txt";
		BufferedWriter bw2 = new BufferedWriter(new FileWriter(out+"2"));
		bw2.flush();bw2.close();
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		String[] modlist={"Acetyl","K+70","K+86","K+84","K+68"};		
		String[] sample={"wolfeiB","wolfeiC","wolfeimC"};
		for (int j=0;j<sample.length;j++){
			for (int i=0;i<modlist.length;i++){
				bw.write(RepPep(paths, modlist[i],sample[j],out+"2"));
			}
		}
		bw.flush();bw.close();
	}
	void process2c() throws IOException{
		String c1="E:/wolfei/20141106/quantArea/20140928wolfeiC.temp";
		String c2="E:/wolfei/20141106/quantArea/20141024wolfeiC.temp";
		String c3="E:/wolfei/20141106/quantArea/20141024wolfeiC_2.temp";
		
		ArrayList<String[]>list=new ArrayList<String[]>();
		String[] a={c1,c2};
		String[] b={c1,c3};
		String[] c={c2,c3};
		list.add(a);list.add(b);list.add(c);
		
		String out="e:/wolfei/20150113-modPeps/ReplicatedPepC.txt";
		BufferedWriter bw2 = new BufferedWriter(new FileWriter(out+"2"));
		bw2.flush();bw2.close();
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		String[] modlist={"Acetyl","K+70","K+86","K+84","K+68"};	
		for (int j=0;j<list.size();j++){
			for (int i=0;i<modlist.length;i++){
				bw.write(RepPep(list.get(j), modlist[i],"wolfeiC",out+"2"));
			}
		}		
		bw.flush();bw.close();
	}
	void process1() throws IOException{
		File folder=new File("E:/wolfei/20141106/quantArea");
		File[] list=folder.listFiles();
		String[] paths=new String[list.length];
		for (int i=0;i<list.length;i++){
			String file=list[i].getAbsolutePath();
			paths[i]=file;
		}
		
		BufferedWriter bw = new BufferedWriter(new FileWriter("e:/wolfei/20150113-modPeps/countShotgun.txt"));
		String[] modlist={"Acetyl,","Butyryl,","3OHbutyryl,","Acetoacetyl,","Crotonyl,"};
		for (int i=0;i<modlist.length;i++){
			bw.write(countSample(paths, modlist[i]));
		}
		
		String[] sample={"wolfeiB","wolfeiC","wolfeimC"};
		for (int j=0;j<sample.length;j++){
			for (int i=0;i<modlist.length;i++){
				bw.write(countRep(paths, modlist[i],sample[j]));
			}
		}
		bw.flush();bw.close();
	}
	private HashSet<String> setExtract(String path,String mod) throws IOException{
		HashSet<String> set=new HashSet<String>();
		BufferedReader br = new BufferedReader(new FileReader(path));
		String line=br.readLine();
		while (line!=null){
			if (line.contains(mod)){
				String[] array=line.split("\t");
				String peptide=array[1]; //extract
				if (peptide.contains(";")){
					String[] p=peptide.split(";");
					for (int i=0;i<p.length;i++){
						set.add(p[i]);
					}
				}else
					set.add(peptide);
			}
			line=br.readLine();
		}
		br.close();
		return set;
	}
	
	//read extract sequences for a biological sample, than compare
	public String countSample(String[] paths, String mod) throws IOException{		
		String returns=mod+"\nsetX\tsetY\tLenX\tLenY\tShareXY\tshareAll\tXYonly\n";
		ArrayList<HashSet<String>> list=new ArrayList<HashSet<String>>();
		String[] sample={"wolfeiB","wolfeiC","wolfeimC"};
		for (int j=0;j<sample.length;j++){
			HashSet<String>aset=new HashSet<String>();
			for (int i=0;i<paths.length;i++){
				if (paths[i].contains(sample[j])&&paths[i].contains("mod2-2")){
					HashSet<String> set=setModPep(paths[i], mod);
					aset.addAll(set);
				}				
			}
			list.add(aset);
		}
		HashSet<String>all=new HashSet<String>(list.get(0)); 
		for (int i=1;i<list.size();i++){
			HashSet<String>temp=new HashSet<String>(all);
			HashSet<String>aset=list.get(i);
			Iterator<String>it=all.iterator();
			while (it.hasNext()){
				String item=it.next();
				if (!aset.contains(item))
					temp.remove(item);
			}
			all.clear();all.addAll(temp);
		}
		int shareall=all.size();
		for (int i=0;i<list.size();i++){
			for (int j=i+1;j<list.size();j++){
				HashSet<String>set1=list.get(i);
				HashSet<String>set2=list.get(j);
				int len1=set1.size(); int len2=set2.size();
				HashSet<String>set12=new HashSet<String>(set1);
				set12.addAll(set2);
				int len=set12.size();
				int common=len1+len2-len, uni1=len1-len, uni2=len2-len;
				returns+="set"+i+"\tset"+j+"\t"+len1+"\t"+len2+"\t"+common+"\t"+shareall+"\t"+(int)(common-shareall)+"\n";
			}
		}			
		return returns;
	}
	
	private HashSet<String> setModPep(String path,String mod) throws IOException{
		HashSet<String> set=new HashSet<String>();
		BufferedReader br = new BufferedReader(new FileReader(path));
		String line=br.readLine();
		while (line!=null){
			if (line.contains(mod)){
				String[] array=line.split("\t");
				String peptide=array[5];
				if (peptide.contains(";")){
					String[] p=peptide.split(";");
					for (int i=0;i<p.length;i++){
						set.add(p[i]);
					}
				}else
					set.add(peptide);
			}
			line=br.readLine();
		}
		br.close();
		return set;
	}
	
	//if ID'd peptides were replicated in different runs
	public String countRep(String[] paths, String mod,String sample) throws IOException{		
		String returns=sample+" "+mod+"\n"+"setX\tsetY\tLenX\tLenY\tShareXY\tshareAll\tXYonly\n";
		ArrayList<HashSet<String>> list=new ArrayList<HashSet<String>>();
		for (int i=0;i<paths.length;i++){
			if (paths[i].contains("mod2-2")&&paths[i].contains(sample)){
				HashSet<String> set=setModPep(paths[i], mod);
				list.add(set);
			}
		}
		HashSet<String>all=new HashSet<String>(list.get(0)); 
		for (int i=1;i<list.size();i++){
			HashSet<String>temp=new HashSet<String>(all);
			HashSet<String>aset=list.get(i);
			Iterator<String>it=all.iterator();
			while (it.hasNext()){
				String item=it.next();
				if (!aset.contains(item))
					temp.remove(item);
			}
			all.clear();all.addAll(temp);
		}
		int shareall=all.size();
		for (int i=0;i<list.size();i++){
			for (int j=i+1;j<list.size();j++){
				HashSet<String>set1=list.get(i);
				HashSet<String>set2=list.get(j);
				int len1=set1.size(); int len2=set2.size();
				HashSet<String>set12=new HashSet<String>(set1);
				set12.addAll(set2);
				int len=set12.size();
				int common=len1+len2-len, uni1=len1-len, uni2=len2-len;
				returns+="set"+i+"\tset"+j+"\t"+len1+"\t"+len2+"\t"+common+"\t"+shareall+"\t"+(int)(common-shareall)+"\n";
			}
		}	
		return returns;
	}
	
	private TreeMap<String,String> ModArea(String path,String mod) throws IOException{
		TreeMap<String, String>tree=new TreeMap<String, String>();
		BufferedReader br = new BufferedReader(new FileReader(path));
		String line=br.readLine(); int dup=0;
		while (line!=null){
			if (line.contains(mod)){
				String[] array=line.split("\t");
				String peptide=array[1];
				String area=array[3];
				String score=array[4];
				if (tree.containsKey(peptide)){
					dup++;
					String[] val=tree.get(peptide).split("\t");
					if (area.length()==0 || val[0].length()==0){
						if (val[0].length()==0)
							tree.put(peptide, area+"\t"+score);
					}else
						System.out.println("rep w/o 0 area len? "+peptide+" "+area +" "+score+" "+tree.get(peptide));
				}else
					tree.put(peptide, area+"\t"+score);
			}
			line=br.readLine();
		}
//		if (dup>0)
//		System.out.println(path+" dup size from temp "+dup);
		br.close();
		return tree;
	}
	
	//if ID'd peptides were replicated in different runs
	public String RepPep(String[] paths, String mod,String sample,String out) throws IOException{	
		BufferedWriter bw=new BufferedWriter(new FileWriter(out,true));
		String returns=sample+" "+mod+"\n"+"peptie\tsetX area\tscore\tsetY area\tscore\tSetZ area\tscore\n";
		ArrayList<TreeMap<String,String>> list=new ArrayList<TreeMap<String,String>>();
		for (int i=0;i<paths.length;i++){
			if (paths[i].contains("temp")&&paths[i].contains(sample)){
				TreeMap<String, String> tree=ModArea(paths[i], mod); //peptide, area \t score
				list.add(tree);
			}
		}
		HashSet<String>all=new HashSet<String>(list.get(0).keySet()); 
		for (int i=1;i<list.size();i++){
			HashSet<String>temp=new HashSet<String>(all);
			HashSet<String>aset=new HashSet<String>(list.get(i).keySet());
			Iterator<String>it=all.iterator();
			while (it.hasNext()){
				String item=it.next();
				if (!aset.contains(item))
					temp.remove(item);
			}
			all.clear();all.addAll(temp);
		}
		bw.write(sample+" "+mod+"\n"+all.size()+"\t");
		for (int i=0;i<list.size();i++){
			bw.write(list.get(i).size()+"\t");
		}
		Iterator<String> it=all.iterator();
		while (it.hasNext()){
			String s=it.next();
			returns+=s+"\t";
			for (int i=0;i<list.size();i++){
				double area=0;
				if (!list.get(i).containsKey(s))
					System.out.println("error not having this key "+s);
				String[] st=list.get(i).get(s).split("\t");
				list.get(i).remove(s);
				if (!st[0].equals(""))
					area=Math.log10(Double.parseDouble(st[0]));
				returns+=area+"\t"+st[1]+"\t";
			}
			returns+="\n";
		}
		returns+="\n";
		for (int i=0;i<list.size();i++){
			TreeMap<String, String> tree=list.get(i);
			bw.write(tree.size()+"\t");
			returns+=sample+"\t"+mod+"\tlist"+i+"\n";
			for (Map.Entry<String, String> en:tree.entrySet()){
				double area=0;
				returns+=en.getKey()+"\t";
				String[] st=en.getValue().split("\t");
				if (!st[0].equals(""))
					area=Math.log10(Double.parseDouble(st[0]));
				returns+=area+"\t"+st[1]+"\n";
			}
			returns+="\n";
		}
		bw.write("\n");
		bw.flush();bw.close();
		return returns;
	}

}
