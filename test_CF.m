% define the device (either Arduino or Seeeduino) serial communication
device = serialport("COM10",9600);

% initialized gait cycle number
    m = 0;
    
while(true)
    % reads the gait cycle number as a char value
    m_char = readline(device);
    % converts to double
    m = str2double(m_char);
   
     if(m>=2)
        % defines new T_d as double, rounded to make sure it is an int
        Td_int = round(m*1000);
        % converts double to string 
        Td = num2str(Td_int);
        % writes string to serial port
        writeline(device,Td);
     end
end
