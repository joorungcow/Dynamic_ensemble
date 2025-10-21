function isResponse = PageTurner(isResponse,kb)
    KbQueueCheck;
    while ~isResponse
        [pressed, firstPress]=KbQueueCheck([]);
        if pressed
            if firstPress(kb.spaceKey)
                isResponse = 1;
            elseif firstPress(kb.escKey)
                %Screen('LoadNormalizedGammaTable', dp.wPtr, repmat(linspace(0,1, 256)',1,3));
                ListenChar(1);
                ShowCursor;
                sca;
                break
            end
        end
    end
    
%     isResponse = 0;
%     while ~isResponse
%         [keyIsDown, secs, keyCode] = KbCheck(-1);
%         if keyCode(kb.spaceKey)
%             isResponse = 1;
%         elseif keyCode(kb.escKey)
%             ListenChar(1);
%             ShowCursor;
%             sca;
%             break
%         end
%     end