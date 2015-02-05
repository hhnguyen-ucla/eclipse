package quant;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

public class quantMod {

	public static void main(String[] arg) throws IOException{
		quantMod test=new quantMod();
		ArrayList<String> path=new ArrayList<>();
		ArrayList<TreeMap<String, String>> list=new ArrayList<>();
		path.add("e:/wolfei/wolfeiB-IP.out");
		path.add("e:/wolfei/wolfeiC-IP.out");
		String out="e:/wolfei/20141010/wolfeiAll-IP.out";
		
//		ArrayList<String> path=new ArrayList<>();
//		ArrayList<TreeMap<String, String>> list=new ArrayList<>();
//		path.add("e:/wolfei/20141001/wolfeiB.temp");
//		path.add("e:/wolfei/20141001/wolfeiC.temp");
//		path.add("e:/wolfei/20141001/wolfeimC.temp");
//		path.add("e:/wolfei/20140911/B.temp");
//		path.add("e:/wolfei/20140911/C.temp");
//		String out="e:/wolfei/20141010/modShotgun.out";
		for (int i=0;i<path.size();i++){
			test.cleanMod(path.get(i),path.get(i)+"2");
			list.add(test.readmod(path.get(i)+"2"));
		}
		test.matchMod(list,out);
	}

	//from .temp2 file: prot pep mod vs area score & area score (non mod)
		public TreeMap<String, String> readmod(String path) throws IOException{
			BufferedReader br = new BufferedReader(new FileReader(path));
			BufferedWriter bw=new BufferedWriter(new FileWriter("e:/wolfei/20141010/test.out",true));
			TreeMap<String, String> tree=new TreeMap<>();
			TreeMap<String, ArrayList<String>> peptree=new TreeMap<>(); //prot pep vs mod
			TreeMap<String, String> normtree=new TreeMap<>();
			String line=br.readLine();
			while (line!=null){
				line=line.replace("\t\t", "\t-\t-");
				String[] array=line.split("\t");
				if (array[2].contains("yl")&&!array[2].contains("Carba")){
					String key=array[0]+"\t"+array[1]+"\t"+array[2]; //prot pep acyl vs area score
					String value=array[3]+"\t"+array[4];					
					tree.put(key, value);
					String pkey=array[0]+"\t"+array[1].replace("*", "").toUpperCase();
					String pval=array[0]+"\t"+array[1]+"\t"+array[2];
					if (peptree.containsKey(pkey)){
						ArrayList<String> list=peptree.get(pkey);
						list.add(pval);
						peptree.put(pkey,list);
					}else{
						ArrayList<String> list=new ArrayList<>();
						list.add(pval);
						peptree.put(pkey,list);
					}
				}else{
					String key=array[0]+"\t"+array[1].replace("c", "C");
					String value=array[3]+"\t"+array[4];
					if (normtree.containsKey(key)){
						String[] varray=normtree.get(key).split("\t");
						double vscore=Math.abs(Double.parseDouble(varray[1]));
						double score=Math.abs(Double.parseDouble(array[4]));
						if (score>vscore)
							normtree.put(key, value);
					}else
						normtree.put(key,value);
				}
				line=br.readLine();
			}
			br.close();			
			
			TreeMap<String, String> returns=new TreeMap<>();
			//now it's time to put nonmod pep together with mod pep
			for (Map.Entry<String, ArrayList<String>> en : peptree.entrySet()){
				ArrayList<String> list=en.getValue();
				String key=en.getKey();
				if (normtree.containsKey(key)){
					for (int i=0;i<list.size();i++){
						String mkey=list.get(i);
						returns.put(mkey, tree.get(mkey)+"\t"+normtree.get(key));
						tree.remove(mkey);
					}
				}
			}
			//for peptides without unmod version
			for (Map.Entry<String, String> en:tree.entrySet()){
				returns.put(en.getKey(), en.getValue()+"\t-\t-");
			}
			
			for (Map.Entry<String, String>e:returns.entrySet()){
				bw.write(e.getKey()+"\t"+e.getValue()+"\n");
			}
			bw.flush();bw.close();
			return returns;
		}

	//if it has non-mod pep, write as well
	public void matchMod(ArrayList<TreeMap<String, String>> list,String out)throws IOException{
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		ArrayList<ArrayList<String>> keylist=new ArrayList<>();
		for (int i=0;i<list.size();i++){
			ArrayList<String> ikeylist=new ArrayList<>();
			for (Map.Entry<String, String> entry:list.get(i).entrySet()){
				ikeylist.add(entry.getKey());
			}
			keylist.add(ikeylist);
		}		

		for (int i=0;i<list.size();i++){
			for (int l=0;l<keylist.get(i).size();l++){
				String key=keylist.get(i).get(l);
				if (!key.equals("")){
					String value=list.get(i).get(key);
					String[] array=value.split("\t");
					int len=array.length;
					//patch the first columns
					bw.write(key+"\t");
					for (int p=0;p<i;p++){
						for (int q=0;q<len;q++){
							bw.write("\t");
						}
					}				
					//check if everyone has it
					for (int j=i;j<list.size();j++){
						TreeMap<String, String> next=list.get(j);
						if (next.containsKey(key)){
							bw.write(next.get(key)+"\t");
							keylist.get(j).remove(key);
							list.get(j).remove(key);
						}else{
							for (int q=0;q<len;q++){
								bw.write("\t");
							}
						}
					}
					bw.write("\n");
				}
			}
		}		
		bw.flush();bw.close();
	}

	//from .temp file
	//test the same ion, c-term is k
	public void cleanMod(String path,String out) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(path));
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		TreeMap<String, String> tree=new TreeMap<>();
		String line=br.readLine();
		while (line!=null){
			line=line.replace("\t\t", "\t-\t-");
			String[] array=line.split("\t");
			if (array[2].contains("K")){ //mod on K
				String p=array[1]; //peptide
				String newp=p;
				boolean proceed=true, innerk=false;
				String l=p.charAt(p.length()-1)+"";
				if (l.equals("k")){
					for (int i=0;i<p.length()-1;i++){ //test last residue
						String c=p.charAt(i)+"";
						if (c.equals("K")||c.equals("k")){
							innerk=true;
							newp=p.substring(0,p.length()-1).replace("K", "k")+"K*";
							break;
						}
					}
					if (innerk==false)
						proceed=false;
				}
				

				if (proceed==true){ //not the c-term residue

					String[] a=array[2].split(","); //mod
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
							acyl.add("3OHButyryl");
						else if (mod.contains("K+84"))
							acyl.add("acetoacetyl");
					}
					Iterator<String> it=acyl.iterator();
					while (it.hasNext())
						acylation+=it.next()+",";

					String key=array[0]+"\t"+newp.toUpperCase().replace("*", "")+"\t"+acylation; //prot pep acyl 
					String value=array[0]+"\t"+newp+"\t"+acylation+"\t"+array[3]+"\t"+array[4];//prot pep acyl area score
					if (tree.containsKey(key)){
						System.out.println(key+" "+value);
						String[] v=tree.get(key).split("\t");
						double vscore=Math.abs(Double.parseDouble(v[4]));
						double score=Math.abs(Double.parseDouble(array[4]));
						if (score>vscore)
							tree.put(key, value);
					}else
						tree.put(key, value);
				}
			}else{
				bw.write(line+"\n");
			}
			line=br.readLine();
		}
		for (Map.Entry<String, String> entry: tree.entrySet()){
			bw.write(entry.getValue()+"\n");
		}
		
		br.close();bw.flush();bw.close();
	}
}

