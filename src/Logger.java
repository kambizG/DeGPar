import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class Logger {

	PrintWriter log;

	public Logger(String logfile) {
		try {
			log = new PrintWriter(new BufferedWriter(new FileWriter(logfile
					+ ".log", true)));
		} catch (IOException ioe) {
			System.err.println(ioe);
			ioe.printStackTrace();
		}
	}

	public String formatEntry(String entry) {
		return new java.util.Date() + " " + entry + "\n";
	}

	public String logEntry(String entry) {
		entry = formatEntry(entry);
		log.println(entry);
		log.flush();
		return entry;
	}
}
