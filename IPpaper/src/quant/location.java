package quant;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class location {
	
	public static void main(String[] arg) throws IOException{
		location.readCello("e:/wolfei/cello.location", "e:/wolfei/cello2.location");
		
	}
	
	public static void readCello(String path,String out) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(path));
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		String line=br.readLine();
		while (line!=null){
//			if (line.contains("*"))
//				System.out.println(line);
			if (line.contains("SeqID")){
				String id=line.replace("SeqID: ","");
				String location="";
				line=br.readLine();
				while (!line.contains("***")){
//					System.out.println(line);
					if (line.contains("*")){
//						System.out.println(line);
						line=line.replace(" \t  ", "");
						location+=line.substring(0,line.indexOf("\t")).trim()+";";
//						location+=line;
					}
					line=br.readLine();
				}
				bw.write(id+"\t"+location+"\n");
			}
			line=br.readLine();
		}
		br.close();bw.flush();bw.close();
	}

}
