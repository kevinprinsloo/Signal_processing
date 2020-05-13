% Import text file
file_path='G:\My Drive\MATLAB\SpeechProduction\transcript\ogg_transcript.txt';
fileID = fopen(file_path,'r');

%% Go through each line
new_section=false;
clc
i=1;
word_count=0;
phrase_count=0;
phrase_end_word=0;
word_confidence=[];
word_in_phrase=[];
speaker_cnt=0;
while ~feof(fileID) 
    tline = fgets(fileID);
    spaces=0;
    j=1;
    while (tline(j)==' ') % check how many spaces before line
        spaces=spaces+1;
        j=j+1;
    end
    if spaces==21 %% word or onset or offset
        word_count=word_count+1;
        word{word_count}=tline(or(or(and(tline>=97,tline<=122),and(tline>=65,tline<=90)),tline==42));
        onset_char_arr = fgets(fileID);
        onset(word_count)=str2num(onset_char_arr);
        offset_char_arr = fgets(fileID);
        offset(word_count)=str2num(offset_char_arr);
        i=i+1; %for onset
        i=i+1; %for offset          
    end
    if spaces==15 % confidence / transcript
        string=tline(or(and(tline>=97,tline<=122),and(tline>=65,tline<=90)));
        if strcmp(string,'confidence')
        	phrase_count=phrase_count+1;
            phrase_confidence(phrase_count)=str2num(tline(spaces+15:end));
            phrase_start_word=phrase_end_word+1;
            phrase_end_word=word_count;
            word_in_phrase(phrase_start_word:phrase_end_word)=phrase_count;        
            word_confidence(phrase_start_word:phrase_end_word)=phrase_confidence(phrase_count);          

        end
        if length(string)>=10
            if strcmp(string(1:10),'transcript')
                phrase_string=tline(spaces+16:end-3);
                phrase{phrase_count}=phrase_string;                
            end
        end
    end
    
    if spaces==9
        if length(tline)>16
            if strcmp((tline(spaces+2:spaces+8)),'speaker')
               % disp(tline)
               speaker_cnt=speaker_cnt+1;
               speaker(speaker_cnt)=str2num(tline(spaces+11:end-3));
            end
        end
    end
    
    i=i+1;
end
fclose(fileID);
%%
transcript_matrix={};
for word_ind=1:word_count
    transcript_matrix{word_ind,1}=word{word_ind};
    transcript_matrix{word_ind,2}=onset(word_ind);
    transcript_matrix{word_ind,3}=offset(word_ind);
    transcript_matrix{word_ind,4}=word_in_phrase(word_ind);
    transcript_matrix{word_ind,5}=word_confidence(word_ind);
    transcript_matrix{word_ind,6}=speaker(word_ind);
    transcript_matrix{word_ind,7}=phrase(word_in_phrase(word_ind));
end

