function showInfo(params,center,dp,kb,iRun)
    if iRun == 0
        %[oldFontName,oldFontNumber,oldTextStyle] = Screen('TextFont', dp.wPtr, '나눔 고딕');
        Screen(dp.wPtr,'TextSize', 30);
        str = double('1) 실험에 참가해주셔서 감사합니다.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy-220);
        str = double('2) 제시되는 자극을 보고 직감적으로 평균 크기가 더 커보이는 쪽의 방향키를 누르십시오.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy-120);
        str = double('3) 연습 시행 30회 후 본 시행을 진행할 것입니다.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy);
        str = double('4) 본 시행 중간에는 아홉 번의 쉬는 시간이 있습니다.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy+120);
        str = double('5) 질문 사항이 없으시면 Space Key를 눌러 계속 진행하십시오.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy+220);
        Screen('Flip', dp.wPtr);   

        isResponse = 0;
        isResponse = PageTurner(isResponse,kb);
            
        make_fixation(dp,center,params.fixSize,[0 0 0]);
        Screen('Flip', dp.wPtr);
        t1 = GetSecs;
        while GetSecs-t1 <= .5; end
        
        str = double('연습 시행에서는 피드백이 제공됩니다.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy-120);
        make_fixation(dp,center,params.fixSize,[0 0 0]);
        str = double('넘기시려면 Space Key를 누르십시오.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy+60);
        Screen('Flip', dp.wPtr); 

        isResponse = 0;
        isResponse = PageTurner(isResponse,kb);
        
        make_fixation(dp,center,params.fixSize,[0 0 0]);
        Screen('Flip', dp.wPtr);
        t1 = GetSecs;
        while GetSecs-t1 <= .5; end
        
        str = double('맞출 경우 아래와 같이 응시점이 초록색으로 바뀝니다.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy-120);
        make_fixation(dp,center,params.fixSize,[0 .8 0]);
        str = double('넘기시려면 Space Key를 누르십시오.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy+60);
        Screen('Flip', dp.wPtr);

        isResponse = 0;
        isResponse = PageTurner(isResponse,kb);
        
        make_fixation(dp,center,params.fixSize,[0 0 0]);
        Screen('Flip', dp.wPtr);
        t1 = GetSecs;
        while GetSecs-t1 <= .5; end
        
        str = double('틀릴 경우 아래와 같이 응시점이 빨간색으로 바뀝니다.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy-120);
        make_fixation(dp,center,params.fixSize,[.8 0 0]);
        str = double('넘기시려면 Space Key를 누르십시오.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy+60);
        Screen('Flip', dp.wPtr);

        isResponse = 0;
        isResponse = PageTurner(isResponse,kb);
        
        make_fixation(dp,center,params.fixSize,[0 0 0]);
        Screen('Flip', dp.wPtr);
        t1 = GetSecs;
        while GetSecs-t1 <= .5; end
        
        str = double('아래의 응시점을 계속 주시해주시길 바랍니다.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy-120);
        make_fixation(dp,center,params.fixSize,params.fixColor)
        str = double('넘기시려면 Space Key를 누르십시오.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy+60);
        Screen('Flip', dp.wPtr);

        isResponse = 0;
        isResponse = PageTurner(isResponse,kb);
        
        make_fixation(dp,center,params.fixSize,[0 0 0]);
        Screen('Flip', dp.wPtr);
        t1 = GetSecs;
        while GetSecs-t1 <= .5; end
        
        str = double('연습 시행을 시작하겠습니다.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy-120);
        make_fixation(dp,center,params.fixSize,params.fixColor)
        str = double('Space Key를 눌러 시작하십시오.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy+60);
        make_fixation(dp,center,params.fixSize,params.fixColor)
        Screen('Flip', dp.wPtr);

    elseif iRun == 1
        str = double('지금부터 본 시행을 시작하겠습니다.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy-220);
        str = double('본 시행에서는 피드백이 제공되지 않습니다.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy-120);
        str = double('Space Key를 누르면 실험이 시작됩니다.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy+60);
        make_fixation(dp,center,params.fixSize,params.fixColor)
        Screen('Flip', dp.wPtr);
    else
        if params.forcedBreak
            t1 = GetSecs;
            while GetSecs-t1 <= 30
                nt = fix(GetSecs-t1);
                %mins = fix((nt)/60);
                %secs = fix(nt-(mins*60));
                str = double('30초 간 휴식해주십시오.');
                strBounds = Screen('TextBounds', dp.wPtr, str);
                Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy-120);
                %str = sprintf('%d분 %d초/3분', mins, secs);
                str = sprintf('%d초 / 30초', nt);
                strBounds = Screen('TextBounds', dp.wPtr, str);
                Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy+60);
                make_fixation(dp,center,params.fixSize,params.fixColor)
                Screen('Flip', dp.wPtr);
                t2 = GetSecs;
                while GetSecs-t2 <= 1; end
            end
        end

        str = sprintf('%d / 10', iRun);
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy-120);
        str = double('Space Key를 누르면 다음 시행이 시작됩니다.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy+60);
        make_fixation(dp,center,params.fixSize,params.fixColor)
        Screen('Flip', dp.wPtr);
    end
    
    isResponse = 0;
    isResponse = PageTurner(isResponse,kb);

end

