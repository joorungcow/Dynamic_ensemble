function showInfo(params,center,dp,kb,iRun)
    if iRun == 0
        %[oldFontName,oldFontNumber,oldTextStyle] = Screen('TextFont', dp.wPtr, '���� ���');
        Screen(dp.wPtr,'TextSize', 30);
        str = double('1) ���迡 �������ּż� �����մϴ�.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy-220);
        str = double('2) ���õǴ� �ڱ��� ���� ���������� ��� ũ�Ⱑ �� Ŀ���̴� ���� ����Ű�� �����ʽÿ�.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy-120);
        str = double('3) ���� ���� 30ȸ �� �� ������ ������ ���Դϴ�.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy);
        str = double('4) �� ���� �߰����� ��ȩ ���� ���� �ð��� �ֽ��ϴ�.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy+120);
        str = double('5) ���� ������ �����ø� Space Key�� ���� ��� �����Ͻʽÿ�.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy+220);
        Screen('Flip', dp.wPtr);   

        isResponse = 0;
        isResponse = PageTurner(isResponse,kb);
            
        make_fixation(dp,center,params.fixSize,[0 0 0]);
        Screen('Flip', dp.wPtr);
        t1 = GetSecs;
        while GetSecs-t1 <= .5; end
        
        str = double('���� ���࿡���� �ǵ���� �����˴ϴ�.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy-120);
        make_fixation(dp,center,params.fixSize,[0 0 0]);
        str = double('�ѱ�÷��� Space Key�� �����ʽÿ�.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy+60);
        Screen('Flip', dp.wPtr); 

        isResponse = 0;
        isResponse = PageTurner(isResponse,kb);
        
        make_fixation(dp,center,params.fixSize,[0 0 0]);
        Screen('Flip', dp.wPtr);
        t1 = GetSecs;
        while GetSecs-t1 <= .5; end
        
        str = double('���� ��� �Ʒ��� ���� �������� �ʷϻ����� �ٲ�ϴ�.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy-120);
        make_fixation(dp,center,params.fixSize,[0 .8 0]);
        str = double('�ѱ�÷��� Space Key�� �����ʽÿ�.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy+60);
        Screen('Flip', dp.wPtr);

        isResponse = 0;
        isResponse = PageTurner(isResponse,kb);
        
        make_fixation(dp,center,params.fixSize,[0 0 0]);
        Screen('Flip', dp.wPtr);
        t1 = GetSecs;
        while GetSecs-t1 <= .5; end
        
        str = double('Ʋ�� ��� �Ʒ��� ���� �������� ���������� �ٲ�ϴ�.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy-120);
        make_fixation(dp,center,params.fixSize,[.8 0 0]);
        str = double('�ѱ�÷��� Space Key�� �����ʽÿ�.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy+60);
        Screen('Flip', dp.wPtr);

        isResponse = 0;
        isResponse = PageTurner(isResponse,kb);
        
        make_fixation(dp,center,params.fixSize,[0 0 0]);
        Screen('Flip', dp.wPtr);
        t1 = GetSecs;
        while GetSecs-t1 <= .5; end
        
        str = double('�Ʒ��� �������� ��� �ֽ����ֽñ� �ٶ��ϴ�.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy-120);
        make_fixation(dp,center,params.fixSize,params.fixColor)
        str = double('�ѱ�÷��� Space Key�� �����ʽÿ�.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy+60);
        Screen('Flip', dp.wPtr);

        isResponse = 0;
        isResponse = PageTurner(isResponse,kb);
        
        make_fixation(dp,center,params.fixSize,[0 0 0]);
        Screen('Flip', dp.wPtr);
        t1 = GetSecs;
        while GetSecs-t1 <= .5; end
        
        str = double('���� ������ �����ϰڽ��ϴ�.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy-120);
        make_fixation(dp,center,params.fixSize,params.fixColor)
        str = double('Space Key�� ���� �����Ͻʽÿ�.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy+60);
        make_fixation(dp,center,params.fixSize,params.fixColor)
        Screen('Flip', dp.wPtr);

    elseif iRun == 1
        str = double('���ݺ��� �� ������ �����ϰڽ��ϴ�.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy-220);
        str = double('�� ���࿡���� �ǵ���� �������� �ʽ��ϴ�.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy-120);
        str = double('Space Key�� ������ ������ ���۵˴ϴ�.');
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
                str = double('30�� �� �޽����ֽʽÿ�.');
                strBounds = Screen('TextBounds', dp.wPtr, str);
                Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy-120);
                %str = sprintf('%d�� %d��/3��', mins, secs);
                str = sprintf('%d�� / 30��', nt);
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
        str = double('Space Key�� ������ ���� ������ ���۵˴ϴ�.');
        strBounds = Screen('TextBounds', dp.wPtr, str);
        Screen('DrawText', dp.wPtr, str, dp.cx-strBounds(3)/2, dp.cy+60);
        make_fixation(dp,center,params.fixSize,params.fixColor)
        Screen('Flip', dp.wPtr);
    end
    
    isResponse = 0;
    isResponse = PageTurner(isResponse,kb);

end

