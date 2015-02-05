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

public class CountMod3 {
	
	public static void main(String[] arg) throws IOException{
		CountMod3 test=new CountMod3();
		String ptt="E:/wolfei/NC_008346.faa";
		File folder=new File("E:/wolfei/20141106-quantArea/quantArea2");
		File[] list=folder.listFiles();
		for (File file:list){
			String path=file.getAbsolutePath();
			if (path.contains(".temp")&&!path.contains("ecoli")&&!path.contains("temp2")){
				System.out.println(path);
				String mod22=path.replace(".temp", ".mod2-2");
				String nomod=path.replace(".temp", ".mod2-nomod");
				String mod23=path.replace(".temp", ".mod2-3");
				test.getUnmodK(path,mod22 ,nomod );
				test.fixupMod22(mod22, nomod, mod23);
			}
		}
	}
	
	private TreeMap<String, String> readnomod(String nomod) throws IOException{
		TreeMap<String, String> tree=new TreeMap<String, String>(); //prot-pepx vs count-score/len:::-sumArea
		BufferedReader br = new BufferedReader(new FileReader(nomod));
		String line=br.readLine();
		while (line!=null){
			String[]array=line.replace("\t\t", "\t-\t-").split("\t");
			String prot=array[0];
			String pepx=array[1];
			String[] scores=array[4].split(":");
			String[] peps=array[5].split(":");
			int count=scores.length;
			int countspec=0;
			String scoreLen="";
			for (int i=0;i<scores.length;i++){
				String s=scores[i];
				String spec=s.substring(s.indexOf("(")+1,s.indexOf(")"));
				countspec+=Integer.parseInt(spec);
				scoreLen+=scores[i]+"|"+peps[i].length()+":";
			}
			scoreLen=scoreLen.substring(0,scoreLen.length()-1); //get rid of :
			String sumArea=array[3].replace("-", "");
			if (sumArea.equals("0.0"))
				sumArea="";
			String key=prot+"\t"+pepx;
			String value="("+countspec+")"+count+"\t"+scoreLen+"\t"+sumArea;
			tree.put(key, value);
			line=br.readLine();
		}
		br.close();
		return tree;
	}
	
	public void fixupMod22(String mod22,String nomod,String out) throws IOException{
		TreeMap<String, String> tree=readnomod(nomod); //prot-pepx vs count-score/len:::-sumArea
		BufferedReader br = new BufferedReader(new FileReader(mod22));
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		String line=br.readLine();
		while (line!=null){
			String[]array=line.replace("\t\t", "\t-\t-").split("\t");
			String prot=array[0];
			String pepx=array[1];
			String acylation=array[2];
			String pos=array[6];
			String[] scores=array[4].split(":");
			String[] peps=array[5].split(":");
			int count=scores.length;
			int countspec=0;
			String scoreLen="";
			for (int i=0;i<scores.length;i++){
				String s=scores[i];
				String spec=s.substring(s.indexOf("(")+1,s.indexOf(")"));
				countspec+=Integer.parseInt(spec);
				scoreLen+=scores[i]+"|"+peps[i].length()+":";
			}
			scoreLen=scoreLen.substring(0,scoreLen.length()-1);
			String sumArea=array[3].replace("-", "");
			String key=prot+"\t"+pepx;
			String[] firstvalue=tree.get(tree.firstKey()).split("\t");
			String value=""; 
			double varea=0.0; //to catch
			String vareastr="";//to parse
			String logvarea="";//to write
			if (tree.containsKey(key)){
				value=tree.get(key);
				String[] varray=value.replace("\t", "\t-").split("\t");
//				System.out.println(value);
				vareastr=varray[2].replace("-", "");
				if (vareastr.length()>0){
					varea=Double.parseDouble(vareastr);
//					logvarea=((double)Math.round(Math.log10(varea)*10)/10)+" log";
				}
			}else{
				for (int i=0;i<firstvalue.length-1;i++)
					value+="\t";
			}
			String towritevalue="";
			
			if (sumArea.equals("0.0")){
				sumArea="";
				String vvalue=value.replace("\t", "");
				if (vvalue.equals("")){
					towritevalue=value+"";
				}else
					towritevalue=value.replace(vareastr, logvarea);
			}else{
				double area=Double.parseDouble(sumArea);
				sumArea=((double)Math.round(Math.log10(area)*10)/10)+"";
				if (varea>0){
					logvarea=((int)Math.round(area/(varea+area)*100))+"";
				}else{ //varea=0.0 or unmod pep doesnt have area
					logvarea="99";
				}
				String vvalue=value.replace("\t", "");
				if (vvalue.equals("")){
//					logvarea="100";
					towritevalue=value+"100";
				}else
					towritevalue=value.replace(vareastr, logvarea);
			}

			bw.write(key+"\t"+pos+"\t"+acylation+"\t("+countspec+")"+count+"\t"+scoreLen+"\t"+sumArea+"\t"+towritevalue+"\n");			
			line=br.readLine();
		}
		br.close();
		bw.flush();bw.close();
		
	}
	
	
	private TreeMap<String, HashSet<String>> makemodtree(String mod22p) throws IOException{
		TreeMap<String, HashSet<String>> modtree=new TreeMap<String, HashSet<String>>(); //prot - pos K
		BufferedReader brm=new BufferedReader(new FileReader(mod22p));
		HashSet<String> modpep=new HashSet<String>();
		String linem=brm.readLine();
		while (linem!=null){
			String[] array=linem.split("\t");
			String prot=array[0];
			String pos=array[6];
			if (modtree.containsKey(prot)){
				HashSet<String>set=modtree.get(prot);
				set.add(pos);
				modtree.put(prot,set);
			}else{
				HashSet<String> set=new HashSet<String>();
				set.add(pos);
				modtree.put(prot, set);
			}
			linem=brm.readLine();
		}
		brm.close();
		return modtree;
	}
	
	//input .temp and mod2-2: read mod2-2 first for hashset of mod pep, then read temp and extract seq from there
	//compare if extract-temp has extract-mod22, write
	public void getUnmodK(String tempp, String mod22p, String out) throws IOException{
		String ptt="E:/wolfei/NC_008346.faa";
		TreeMap<String, String>prottree=protseq(ptt, modprot(tempp));
		TreeMap<String, HashSet<String>> modtree=makemodtree(mod22p); //prot - pos K
		TreeMap<String, String> tree=new TreeMap<String, String>();

		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		BufferedReader brt=new BufferedReader(new FileReader(tempp));
		String linet=brt.readLine();
		while (linet!=null){
			linet=linet.replace("\t\t", "\t-\t-");
			String[] array=linet.split("\t");
			String prot=array[0].replace("-", "");
			if (modtree.containsKey(prot)){			
				String pep=array[1].replace("-", "");
				String mods=array[2].replace("-", "");
				String area=array[3].replace("-", "");
				String scores=array[4].replace("-", "");
				String numspec=array[5].replace("-", "");
				if (scores.contains("."))
					scores=scores.substring(0,scores.indexOf("."));
				scores="("+numspec+")"+scores;
				HashSet<String>modpos=modtree.get(prot);
				if (pep.contains("K")&&!mods.contains("K")){ //non-mod pep containing K
					String key="",value="";
					double darea=0.0;
					if (area.length()>0)
						darea=Double.parseDouble(area);
					if (prot.contains(";")){
						String[] protarray=prot.split(";");
						for (int i=0;i<protarray.length;i++){
							String protseq=prottree.get(protarray[i]);
							ArrayList<String>Kpos=extractKPos(pep, protseq);
							for (int j=0;j<Kpos.size();j++){
								String[] a=Kpos.get(j).split("\t");
								String pos=a[0], extract=a[1];
								if (modpos.contains(pos)){									
									key=prot+"\t"+extract+"\t"; //prot extractpep acyl vs area score pep
									value=darea+"\t"+scores+"*\t"+pep+"*\t"+pos;
									if (tree.containsKey(key)){
										String[] v=tree.get(key).split("\t");
										String vscore=v[1];
										double varea=Double.parseDouble(v[0]);
										varea=varea+darea;
										String vpep=v[2];
										value=varea+"\t"+vscore+":"+scores+"*\t"+vpep+":"+pep+"*\t"+pos;
										tree.put(key, value);
									}else
										tree.put(key, value);
								}
							}
						}
					}else{
						String protseq=prottree.get(prot);
						ArrayList<String>Kpos=extractKPos(pep, protseq);
						for (int j=0;j<Kpos.size();j++){
							String[] a=Kpos.get(j).split("\t");
							String pos=a[0], extract=a[1];
							if (modpos.contains(pos)){
								key=prot+"\t"+extract+"\t"; //prot extractpep acyl vs area score pep
								value=darea+"\t"+scores+"\t"+pep+"\t"+pos;
								if (tree.containsKey(key)){
									String[] v=tree.get(key).split("\t");
									String vscore=v[1];
									double varea=Double.parseDouble(v[0]);
									varea=varea+darea;
									String vpep=v[2];
									value=varea+"\t"+vscore+":"+scores+"\t"+vpep+":"+pep+"\t"+pos;
									tree.put(key, value);
								}else
									tree.put(key, value);
							}
						}
					}
				}	
			}
			linet=brt.readLine();
		}
		for (Map.Entry<String, String> e:tree.entrySet()){
			bw.write(e.getKey()+"\t"+e.getValue()+"\n");
		}
		bw.flush();bw.close();brt.close();
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
		return protset;
	}
	
	private ArrayList<String> extractKPos(String pep,String seq){
		pep=pep.toUpperCase();
		ArrayList<String>list=new ArrayList<String>();
		int peppos=seq.indexOf(pep);
		for (int i=0;i<pep.length();i++){
			String c=pep.charAt(i)+"";
			if (c.equals("K")){
				int Kpos=peppos+i;
				int start=Kpos-5;
				int end=Kpos+6;
//				System.out.println(peppos+" "+start+" "+end+" "+seq.length()+" "+pep);
				if (start<0)
					start=0;
				if (end>=seq.length())
					end=seq.length();
				String extract=seq.substring(start,end);
				list.add(Kpos+"\t"+extract);
			}
		}
		return list;
	}
}
