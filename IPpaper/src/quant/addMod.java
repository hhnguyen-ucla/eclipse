package quant;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map;
import java.util.TreeMap;

public class addMod {
	
	public static void main(String[] arg) throws IOException{
		addMod test=new addMod();
		File folder=new File("E:/wolfei/20141106/quantArea");
		File[] list=folder.listFiles();
		for (File file:list){
			String path=file.getAbsolutePath();
			if (path.contains("mod2-2")){
				test.redisplayMod(path, path.replace("2-2", "4"));
			}
		}
		
//		String path="e:/wolfei/20141106/mod.all";
//		String ptitle="e:/wolfei/20141106/mod-title.all";
//		test.modptt(path, ptitle, path.replace("all","ptt"));
	}
	
	
	private TreeMap<String, String> readtitle(String ptitle) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(ptitle));
		TreeMap<String, String> tree=new TreeMap<>();
		String line=br.readLine();
		while (line!=null){
			String[] array=line.split("\t");
			String pid=array[0];
			String value=line.replace(array[0]+"\t","");
			tree.put(pid, value);
			line=br.readLine();
		}
		br.close();
		return tree;
	}
	public void modptt(String path,String ptitle,String out) throws IOException{
		TreeMap<String, String> tree=readtitle(ptitle);
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		BufferedReader br = new BufferedReader(new FileReader(path));
		String line=br.readLine();
		while (line!=null){
			String[] array=line.split("\t");
			String pid=array[0];
			for (int i=4;i<array.length;i+=2){
				String area=array[i];
				if (area.length()>6&&!area.contains(" ")){
					double logarea=(double)Math.round(Math.log10(Double.parseDouble(area))*10)/10;
					line=line.replace(area, logarea+"");
				}
				
			}
			
			bw.write(line+"\t"+tree.get(pid)+"\n");
			line=br.readLine();
		}
		br.close();bw.flush();bw.close();
	}
	
	public void redisplayMod(String path,String out) throws IOException{
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		BufferedReader br = new BufferedReader(new FileReader(path));
		String line=br.readLine();
		TreeMap<String, TreeMap<String, String>> tree=new TreeMap<>(); //pid-pep-pos vs (amod vs score-area)
		while (line!=null){
			String[] array=line.split("\t");
			if (!array[1].equals("")&&!array[0].equals("")){
				String pos=array[6];
				if (pos.contains(","))
					pos=pos.substring(0,pos.indexOf(","));
				String key=array[0]+","+array[1]+","+pos;
				String mkey=array[2];
				if (array[3].equals("0.0"))
					array[3]="";
				array[4]=array[4].replace(":", ";");
				String mvalue=array[4]+"\t"+array[3];
				if (tree.containsKey(key)){				
					TreeMap<String, String> mod=tree.get(key);
					mod.put(mkey, mvalue);
					tree.put(key, mod);
				}else{
					TreeMap<String, String> mod=new TreeMap<>();
					mod.put(mkey, mvalue);
					tree.put(key, mod);
				}
			}
			
			line=br.readLine();
		}
		
		for (Map.Entry<String, TreeMap<String, String>> e:tree.entrySet()){
			String key=e.getKey();
			TreeMap<String, String> mod=e.getValue();
			bw.write(key+"\t");
			String[] modlist={"Acetyl,","Butyryl,","3OHbutyryl,","Acetoacetyl,","Crotonyl,"};
			for (int i=0;i<modlist.length;i++){
				if (mod.containsKey(modlist[i])){
					bw.write(mod.get(modlist[i])+"\t");
				}else
					bw.write("\t\t");
			}
			bw.write("\n");
		}
		br.close();bw.flush();bw.close();
	}

}
