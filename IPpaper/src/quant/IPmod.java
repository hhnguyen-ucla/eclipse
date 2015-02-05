package quant;

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

public class IPmod {
	
	public static void main(String[] arg) throws IOException{
		IPmod test=new IPmod();
		String ptt="E:/wolfei/NC_008346.faa";
		File folder=new File("E:/wolfei/20150113-modPeps/quantArea");
		File[] list=folder.listFiles();
		for (File file:list){
			String path=file.getAbsolutePath();
			if (path.contains(".temp")&&!path.contains("ecoli")){
				System.out.println(path);
//				test.readmod(path, ptt, path.replace("temp", "mod2-2"));
//				test.count(path.replace("temp", "mod2-2"), path.replace("temp", "mod3"));
				test.convertArea(path.replace("temp", "mod2-2"), path.replace("temp","mod2-22"));
			}
		}
	}
	
	private double smallestArea(String path) throws IOException{
		double min=1000000000000.0;String id="";
		BufferedReader br = new BufferedReader(new FileReader(path));
		String line=br.readLine();
		while (line!=null){
			String[] array=line.split("\t");
			double area=Double.parseDouble(array[3]);
			if (area>1&&area<min){
				min=area;
				id=line;
			}
			line=br.readLine();
		}
		br.close();
//		System.out.println(min+" "+id);
		return min;
	}
	//convert areas in mod2-2 into ratios to the lowest one
	public void convertArea(String path,String out) throws IOException{
		double min=smallestArea(path);
		BufferedWriter bw =  new BufferedWriter(new FileWriter(out));
		BufferedReader br = new BufferedReader(new FileReader(path));
		String line=br.readLine();
		while (line!=null){
			String[] array=line.split("\t");
			double area=Double.parseDouble(array[3]);
			int ratio=(int)(Math.round(area/min));
			bw.write(line+"\t"+ratio+"\n");
			line=br.readLine();
		}
		br.close();bw.flush();bw.close();
	}
	
	
	//input .faa file and protset
	private TreeMap<String, String>protseq(String path,HashSet<String>protset) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(path));
		String line=br.readLine();
		TreeMap<String, String> prottree=new TreeMap<>(); //pid, seq
		while (line!=null){
			if (line.contains(">")){
				String[] array=line.replace("|", ";").split(";");
				String pid=array[1],seq="";
				line=br.readLine();
				if (protset.contains(pid)){
					while (line!=null&&!line.contains(">")){
						seq+=line;
						line=br.readLine();
					}
					prottree.put(pid,seq);
				}
			}
			if (line!=null && !line.contains(">"))
				line=br.readLine();
		}
		br.close();
		return prottree;
	}
	
	
	//input .temp file
	private HashSet<String> modprot(String path) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(path));
		String line=br.readLine();
		HashSet<String> protset=new HashSet<>();
		while (line!=null){
			String[] array=line.split("\t");
//			if (array[1].contains("k")&&!array[0].contains("P")&&!array[0].contains("Q")&&!array[0].contains("spike")){
				if (array[1].contains("k")&&array[0].contains("11456")){
				if (array[0].contains(";")){
					String[] a=array[0].split(";");
					for (int i=0;i<a.length;i++)
						protset.add(a[i]);
				}else
					protset.add(array[0]);
			}
			line=br.readLine();
		}
		br.close();
//		System.out.println(protset.size());
		return protset;
	}
	
	//only work for pep with 1 mod
	private String extractPep(String pep,String seq){
		String extractpep="";
		int locak=0,locapep=0;
		for (int i=0;i<pep.length();i++){
			String c=pep.charAt(i)+"";
			if (c.equals("k")){
				locak=i;
				break;
			}
		}
		String pepUpper=pep.toUpperCase();
//		System.out.println(pep+"\n"+seq);
		locapep=seq.indexOf(pepUpper);
		int res=locapep+locak;
		int start=res-5,end=res+6;
		if (start<0)
			start=0;
		if (end>=seq.length())
			end=seq.length();
		extractpep=seq.substring(start,end)+"\t"+res;
		return extractpep;
	}
	
	//from .temp file
	public void readmod(String path,String ptt,String out) throws IOException{
		TreeMap<String, String>prottree=protseq(ptt, modprot(path));
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		BufferedReader br = new BufferedReader(new FileReader(path));
		TreeMap<String, String> tree=new TreeMap<>();
		String line=br.readLine();
		while (line!=null){
			line=line.replace("\t\t", "\t-\t-");
//			System.out.println(line);
			String[] array=line.split("\t");
			String prot=array[0].replace("-", "");
			String pep=array[1].replace("-", "");
			String mods=array[2].replace("-", "");
			String area=array[3].replace("-", "");
			String scores=array[4].replace("-", "");
			if (scores.contains("."))
				scores=scores.substring(0,scores.indexOf("."));
			if (mods.contains("K")&&prot.contains("11456")){
				String extractpep="";
				if (prot.contains(";")){
					String[] protarray=prot.split(";");
					for (int i=0;i<protarray.length;i++){
						String protseq=prottree.get(protarray[i]);
						extractpep+=extractPep(pep, protseq)+"\t";	
					}
				}else
					extractpep=extractPep(pep, prottree.get(prot));			
				String[] a=mods.split(",");
				String acylation="";
				HashSet<String> acyl=new HashSet<>();
				for (int i=0;i<a.length;i++){
					String mod=a[i];
					if (mod.contains("Acetyl"))
						acyl.add("Acetyl");
					else if (mod.contains("K+70"))
						acyl.add("Butyryl");
					else if (mod.contains("K+68"))
						acyl.add("Crotonyl");
					else if (mod.contains("K+86"))
						acyl.add("3OHbutyryl");
					else if (mod.contains("K+84"))
						acyl.add("Acetoacetyl");
				}
				Iterator<String> it=acyl.iterator();
				while (it.hasNext())
					acylation+=it.next()+",";
				
				//not account for more than 1 mod per peptide
				String[] extractarray=extractpep.split("\t");
				
				String key="",value="";
				double darea=0.0;
				if (area.length()>0)
					darea=Double.parseDouble(area);
				if (prot.contains(";")){ //split into 2 entries, add * to score and ID'd peptide
					String[] protarr=prot.split(";");
					for (int i=0;i<extractarray.length;i+=2){
						key=protarr[i/2]+"\t"+extractarray[i]+"\t"+acylation; //prot extractpep acyl vs area score pep
						value=darea+"\t"+scores+"*\t"+pep+"*\t"+extractarray[i+1];
						if (tree.containsKey(key)){
							String[] v=tree.get(key).split("\t");
							String vscore=v[1];
							double varea=Double.parseDouble(v[0]);
							varea=varea+darea;
							String vpep=v[2];
							value=varea+"\t"+vscore+":"+scores+"*\t"+vpep+":"+pep+"*\t"+extractarray[i+1];
							tree.put(key, value);
						}else
							tree.put(key, value);
					}
				}else{				
					key=prot+"\t"+extractarray[0]+"\t"+acylation; //prot extractpep acyl vs area score pep
					value=darea+"\t"+scores+"\t"+pep+"\t"+extractarray[1];
					if (tree.containsKey(key)){
						String[] v=tree.get(key).split("\t");
						String vscore=v[1];
						double varea=Double.parseDouble(v[0]);
						varea=varea+darea;
						String vpep=v[2];
						value=varea+"\t"+vscore+":"+scores+"\t"+vpep+":"+pep+"\t"+extractarray[1];
						tree.put(key, value);
					}else
						tree.put(key, value);
				}
				
			}
			line=br.readLine();
		}
		br.close();
		for (Map.Entry<String, String> e:tree.entrySet()){
			bw.write(e.getKey()+"\t"+e.getValue()+"\n");
		}
		bw.flush();bw.close();
	}
	
	
	public void count(String path,String out) throws IOException{
		TreeMap<String, ArrayList<Integer>> tree2=readMod22(path);
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		TreeMap<String, TreeMap<String, ArrayList<Integer>>> tree=new TreeMap<String, TreeMap<String,ArrayList<Integer>>>();
		for (Map.Entry<String, ArrayList<Integer>> e:tree2.entrySet()){
			String[] array=e.getKey().split("\t");
			String prot=array[0],mod=array[1];
			ArrayList<Integer>counts=e.getValue();

			if (tree.containsKey(prot)){
				TreeMap<String, ArrayList<Integer>> atree=tree.get(prot);
				atree.put(mod,counts);
				tree.put(prot, atree);
			}else{
				TreeMap<String, ArrayList<Integer>> atree=new TreeMap<String, ArrayList<Integer>>();
				atree.put(mod,counts);
				tree.put(prot, atree);
			}
		}
		
		String[] modlist={"Acetyl,","Butyryl,","3OHbutyryl,","Acetoacetyl,","Crotonyl,"};

		for (Map.Entry<String,TreeMap<String, ArrayList<Integer>>> e:tree.entrySet()){
			String prot=e.getKey();
			String towrite=prot+"\t";
			TreeMap<String, ArrayList<Integer>>atree=e.getValue();
			for (int i=0;i<modlist.length;i++){
				if (atree.containsKey(modlist[i])){
					ArrayList<Integer>counts=atree.get(modlist[i]);
					String count=counts.get(0)+"";
					if (counts.get(1)>0)
						count+=";"+counts.get(1);
					towrite+=(count+"\t");
				}else{
					towrite+="\t";
				}
			}
			bw.write(towrite.substring(0,towrite.length()-1)+"\n");
		}
		bw.flush();bw.close();
	}
	
	//need to work on temp file, convert everything to 11 residues
	public TreeMap<String, ArrayList<Integer>> readMod22(String path) throws IOException{
		TreeMap<String, ArrayList<Integer>> tree2=new TreeMap<>(); //prot-mod vs count(unique, share)
		BufferedReader br = new BufferedReader(new FileReader(path));
		String line=br.readLine();
		while (line!=null){
			String[] array=line.replace("\t\t", "\t-\t-").split("\t");
			String prot=array[0].replace("-", "");
			String mod=array[2].replace("-", "");
			String score=array[4];
			String key=prot+"\t"+mod;
			if (tree2.containsKey(key)){
				ArrayList<Integer> alist=tree2.get(key);
				if (score.contains("*")){
					int share=alist.get(1);
					alist.set(1, share+1);
				}else{
					int unique=alist.get(0);
					alist.set(0, unique+1);
				}					
				tree2.put(key, alist);
			}else{
				ArrayList<Integer> alist=new ArrayList<>();
				alist.add(0); alist.add(0);
				if (score.contains("*")){
					int share=alist.get(1);
					alist.set(1, share+1);
				}else{
					int unique=alist.get(0);
					alist.set(0, unique+1);
				}					
				tree2.put(key, alist);
			}
			
			line=br.readLine();
		}
		br.close();
		return tree2;
	}

}
