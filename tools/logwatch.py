import time
import smtplib
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler

import os




#LOG_FILE = "test.log"

LOG_FILE_path = "/hdd/_broker_backend_with_ANTARES/backend_logfiles/start_brokerbackend.log"
LOG_FILE = "start_brokerbackend.log"

#LOG_FILE_path = "/home/nhernits/Documents/_uni/_Antofagasta/_my_science/_LSST/_POI_Variable_Alerts/_conntect_ANTARES/_logfilewatcher/test.log"
#LOG_FILE = "test.log"


class LogHandler(FileSystemEventHandler):
    def __init__(self):
        super().__init__()
        #self.filename = os.path.basename(LOG_FILE_path)
        # Start at the end of the file so we skip old entries
        with open(LOG_FILE_path, "r") as f:
            #print('opened')
            f.seek(0, 2)  # move to end
            self._position = f.tell()
            

    def on_modified(self, event):
        #print('on_modified')
        #print(event.src_path)
        

        # filename = os.path.basename(LOG_FILE_path)
        #print(self.filename)    
        
        if event.src_path.endswith(LOG_FILE): #<<<<<-
            #print('ends with')

            with open(LOG_FILE_path, "r") as f:
                f.seek(self._position)  # jump to last read position
                new_lines = f.readlines()
                #print(new_lines)
                self._position = f.tell()  # update position
                print(new_lines)
                for line in new_lines:
                    if "error " in line.lower() or "https" in line.lower() or "timed out" in line.lower(): # adjust as needed
                 #       print(f"Error detected: {line.strip()}")
                        send_email(line.strip())

def send_email(message):
    #print("Error detected: Sending Email...")

    #print(mailaddress)
    #print(mailpwd)
    
    try:
        with smtplib.SMTP_SSL('smtp.gmail.com', 465) as server:
            server.login(mailaddress, mailpwd)  # Use App Password
            server.sendmail(
            mailaddress, 
            "mail@mail.de", 
            f"Subject: Log Alert\n\n{message}",
            )
        #print('e-mail was sent')

    except smtplib.SMTPException as e:
        print(f"Error: {e}")
        
        
def main():
    
    global mailaddress
    global mailpwd

    keyFile = open('keys.txt', 'r')
    mailaddress = keyFile.readline().rstrip()
    mailpwd = keyFile.readline().rstrip()
    
    keyFile.close()


    print(mailaddress)
    print(mailpwd)
      
    observer = Observer()
    observer.schedule(LogHandler(), LOG_FILE_path, recursive=False)
    observer.start()

    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        observer.stop()
    observer.join()





if __name__ == "__main__":
    main()  
  
  
  

